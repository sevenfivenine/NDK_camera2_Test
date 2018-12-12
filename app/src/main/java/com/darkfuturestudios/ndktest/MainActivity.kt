package com.darkfuturestudios.ndktest

import android.content.Context
import android.content.pm.PackageManager
import android.graphics.Bitmap
import android.graphics.SurfaceTexture
import android.hardware.Camera
import android.hardware.camera2.*
import android.os.*
import android.support.v4.app.ActivityCompat
import android.support.v7.app.AppCompatActivity
import android.util.Log
import android.util.Range
import android.view.TextureView
import android.view.View
import android.widget.SeekBar
import android.widget.Toast
import com.darkfuturestudios.ndktest.CameraController.*
import kotlinx.android.synthetic.main.activity_camera.*
import kotlinx.android.synthetic.main.activity_main.*
import java.io.File
import java.io.FileOutputStream
import java.io.IOException
import java.text.SimpleDateFormat
import java.util.*
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.min
import kotlin.math.round

class MainActivity : AppCompatActivity() {

    //region companion

    companion object {

        const val TAG: String = "MainActivity"
        const val PERMISSION_REQUEST_CAMERA: Int = 7000

        const val SEEK_BAR_EXPOSURE = 0
        const val SEEK_BAR_FOCUS = 1
        const val SEEK_BAR_GAIN = 2
        const val SEEK_BAR_RES = 3

        const val STACK_EXPOSURE_SECONDS = 10.0

        // Used to load the 'native-lib' library on application startup.
        init {
            System.loadLibrary("native-lib")
        }

    }

    //endregion

    //region members

    private lateinit var textureView: TextureView
    private lateinit var fabTakePhoto: View
    private lateinit var seekBars: HashMap<Int, SeekBar>
    private lateinit var surfaceTextureListener: TextureView.SurfaceTextureListener
    private lateinit var cameraManager: CameraManager
    private lateinit var cameraId: String
    private lateinit var galleryFolder: File

    private var cameraController: CameraController? = null
    private var hardwareSupportsCamera2: Boolean = true

    // Camera settings
    private var exposure: Long = 1000000L // camera2
    private var exposureCompensation: Double = 0.0 // camera
    private var focus: Float = 0.0f
    private var gain: Int = 0 // camera2
    private var gainString: String = "" // camera
    private var resolution: CameraController.Size? = null

    //endregion

    //region lifecycle

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_camera)

        // Find views
        textureView = findViewById(R.id.texture_view)
        fabTakePhoto = findViewById(R.id.fab_take_photo)

        // Create seekBars hash map
        seekBars = hashMapOf(SEEK_BAR_EXPOSURE to seek_bar_exposure,
                SEEK_BAR_FOCUS to seek_bar_focus,
                SEEK_BAR_GAIN to seek_bar_gain,
                SEEK_BAR_RES to seek_bar_res)

        for ((key, seekBar) in seekBars) {
            // Set listener
            seekBar.setOnSeekBarChangeListener(object : SeekBar.OnSeekBarChangeListener {
                override fun onProgressChanged(pSeekBar: SeekBar?, progress: Int, fromUser: Boolean) {
                    calculateCameraSetting(progress, key)
                }

                override fun onStartTrackingTouch(pSeekBar: SeekBar) {
                }

                override fun onStopTrackingTouch(pSeekBar: SeekBar) {
                }
            })
        }

        fabTakePhoto.setOnClickListener {
            takePhoto()
        }

        fab_stack.setOnClickListener {
            stackExposure(STACK_EXPOSURE_SECONDS)
        }

        switch_hide_preview.setOnCheckedChangeListener { compoundButton, isChecked ->
            if (isChecked) {
                textureView.visibility = View.INVISIBLE
                image_view_stack.visibility = View.VISIBLE
            } else {
                textureView.visibility = View.VISIBLE
                image_view_stack.visibility = View.INVISIBLE
            }
        }

        cameraManager = applicationContext.getSystemService(Context.CAMERA_SERVICE) as CameraManager

        surfaceTextureListener = object : TextureView.SurfaceTextureListener {
            override fun onSurfaceTextureAvailable(surfaceTexture: SurfaceTexture, width: Int, height: Int) {
                Log.d(TAG, "SurfaceTextureListener: AVAILABLE")
                initializeCamera()
            }

            override fun onSurfaceTextureSizeChanged(surfaceTexture: SurfaceTexture, width: Int, height: Int) {
                //Log.d(TAG, "SurfaceTextureListener: SIZE CHANGED")
            }

            override fun onSurfaceTextureDestroyed(surfaceTexture: SurfaceTexture): Boolean {
                //Log.d(TAG, "SurfaceTextureListener: DESTROYED")
                return false
            }

            override fun onSurfaceTextureUpdated(surfaceTexture: SurfaceTexture) {
                //Log.d(TAG, "SurfaceTextureListener: UPDATED")
            }
        }
    }

    override fun onResume() {
        super.onResume()
        if (textureView.isAvailable) {
            Log.d(TAG, "Texture view is available")
            initializeCamera()
        } else {
            textureView.surfaceTextureListener = surfaceTextureListener
            Log.d(TAG, "Texture view is NOT available")
        }
    }

    override fun onStop() {
        super.onStop()
        closeCamera()
    }

    //endregion

    //region camera

    /**
     * Sets up camera, creates CameraController, starts preview, then calculates initial camera
     * settings
     */
    private fun initializeCamera() {
        val surfaceTexture = textureView.surfaceTexture

        // User has already given access (if not, permissions dialog callback routes back here)
        if (checkPermissions()) {
            setUpCamera()
            createImageGallery()
            cameraController = createCameraController()
            cameraController?.setPreviewTexture(surfaceTexture)
            NDKTestUtil.setCameraDisplayOrientation(this, cameraController, hardwareSupportsCamera2)

            // These will be overwritten once calculateCameraSetting() is called below
            // TODO these values may cause crashes on some devices
            cameraController?.setPreviewSize(1920, 1080)
            cameraController?.setPictureSize(1080, 1920)

            cameraController?.startPreview()

            for ((key, seekBar) in seekBars) {
                // Initialize the camera settings
                calculateCameraSetting(seekBar.progress, key)
            }
        }
    }

    /**
     * Finds the correct cameraId to use (back-facing, not wide angle) and determines if the
     * software and hardware is capable of supporting camera2
     */
    private fun setUpCamera() {
        try {
            /** In case there is more than one back-facing camera, the one with larger focal lengths
             * is the standard camera, while the other is the wide angle lens. Therefore we need to
             * compare focal lengths when iterating through cameras
             *
             * Focal lengths are FloatArrays which often have a length of 1
             * In case they have multiple elements, it makes sense to compare the max from each */
            var maxFocalLengths: FloatArray = floatArrayOf(0.0f)

            val camera2Manager = CameraControllerManager2(this)

            for (cameraId in cameraManager.cameraIdList) {
                val cameraCharacteristics = cameraManager.getCameraCharacteristics(cameraId)
                val focalLengths = cameraCharacteristics.get(CameraCharacteristics.LENS_INFO_AVAILABLE_FOCAL_LENGTHS)
                Log.d(TAG, "Focal lengths available: $focalLengths")
                if (cameraCharacteristics.get(CameraCharacteristics.LENS_FACING) == CameraCharacteristics.LENS_FACING_BACK
                        && focalLengths.max()!! > maxFocalLengths.max()!!) {
                    this.cameraId = cameraId
                    maxFocalLengths = focalLengths
                }

                if (!camera2Manager.allowCamera2Support(cameraId.toInt()))
                    hardwareSupportsCamera2 = false
            }

            Log.d(TAG, "Camera: $cameraId")
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

    /**
     * Creates a gallery in the phones Pictures directory for this app, if not already present
     */
    private fun createImageGallery() {
        // On any Android phone, there is a public directory for pictures
        val storageDirectory = Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_PICTURES)
        galleryFolder = File(storageDirectory, resources.getString(R.string.app_name))
        if (!galleryFolder.exists()) {
            val wasCreated = galleryFolder.mkdirs()
            if (!wasCreated) {
                Log.e(TAG, "Failed to create directories")
            }
        }
    }

    /**
     * Creates and returns the image file captured, to be saved to storage
     */
    @Throws(IOException::class)
    private fun createImageFile(): File {
        val timeStamp = SimpleDateFormat("yyyyMMdd_HHmmss", Locale.getDefault()).format(Date())
        val imageFileName = "image_" + timeStamp + "_"
        return File.createTempFile(imageFileName, ".jpg", galleryFolder)
    }

    /**
     * Create a CameraController
     * Depending on the software (Android 5.0+?) and hardware (is camera2 supported?) this will
     * be either CameraController1 or CameraController2
     */
    private fun createCameraController(): CameraController? {
        var cameraControllerLocal: CameraController?

        try {
            val cameraErrorCallback = CameraController.ErrorCallback {
                if (cameraController != null) {
                    cameraController = null
                }
            }

            val useCamera2 = Build.VERSION.SDK_INT >= Build.VERSION_CODES.LOLLIPOP && hardwareSupportsCamera2

            cameraControllerLocal = if (useCamera2) {
                val previewErrorCallback = CameraController.ErrorCallback {}
                CameraController2(this, cameraId.toInt(), previewErrorCallback, cameraErrorCallback)
            } else {
                CameraController1(cameraId.toInt(), cameraErrorCallback)
            }
        } catch (e: CameraControllerException) {
            e.printStackTrace()
            cameraControllerLocal = null
        }

        return cameraControllerLocal
    }

    /**
     * Closes the CameraController camera
     */
    private fun closeCamera() {
        val cameraControllerLocal = cameraController
        if (cameraController != null) {
            Log.d(TAG, "Closing camera")
            cameraController = null
            cameraControllerLocal?.stopPreview()
            cameraControllerLocal?.release()
        }
    }

    /**
     * Checks for permissions (camera, write external storage)
     * @return True if permissions have already been granted, false if not (will prompt user)
     */
    private fun checkPermissions(): Boolean {
        try {
            Log.d(TAG, "Checking Permissions for Camera use...")
            if (ActivityCompat.checkSelfPermission(this, android.Manifest.permission.CAMERA)
                    == PackageManager.PERMISSION_GRANTED && ActivityCompat.checkSelfPermission(
                            this, android.Manifest.permission.WRITE_EXTERNAL_STORAGE)
                    == PackageManager.PERMISSION_GRANTED) {
                return true
            } else {
                Log.d(TAG, "Permissions Denied!")
                ActivityCompat.requestPermissions(this@MainActivity,
                        arrayOf(android.Manifest.permission.CAMERA,
                                android.Manifest.permission.WRITE_EXTERNAL_STORAGE),
                                PERMISSION_REQUEST_CAMERA)
            }
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

        return false
    }

    /**
     * Take Photo by calling CameraController.takePicture()
     */
    private fun takePhoto() {
        Log.d(TAG, "Take photo")

        if (!textureView.isAvailable) {
            Log.d(TAG, "Texture view not yet available")
            return
        }

        val pictureCallback = object : CameraController.PictureCallback {
            override fun onStarted() {
                Log.d(TAG, "PictureCallback.onStarted()")
            }

            override fun onCompleted() {
                Log.d(TAG, "PictureCallback.onCompleted()")

                /** For camera2, we need to pause the preview to indicate to the user a photo has
                 *  been captured. For camera, pause is done automatically */

                // Pause preview
                if (cameraController is CameraController2) {
                    cameraController?.stopPreview()
                }

                // Start preview
                try {
                    cameraController?.startPreview()
                } catch (e: CameraControllerException) {
                    e.printStackTrace()
                }

            }

            override fun onPictureTaken(data: ByteArray?) {
                Log.d(TAG, "PictureCallback.onPictureTaken()")

                var outputPhoto: FileOutputStream? = null
                try {
                    outputPhoto = FileOutputStream(createImageFile())
                    /** OpenCamera uses a much more sophisticated method of saving images
                     *  This is much simpler, but less versatile
                     *  We don't even use the data argument, just capture from textureView instead
                     */
                    textureView.bitmap.compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
                } catch (e: Exception) {
                    e.printStackTrace()
                } finally {
                    try {
                        outputPhoto?.close()
                    } catch (e: IOException) {
                        e.printStackTrace()
                    }
                }
            }

            override fun onRawPictureTaken(raw_image: RawImage?) {}

            override fun onBurstPictureTaken(images: MutableList<ByteArray>?) {}

            override fun onFrontScreenTurnOn() {}

        }

        val errorCallback = CameraController.ErrorCallback { Log.e(TAG, "Error from takePicture()") }

        cameraController?.takePicture(pictureCallback, errorCallback)
    }

    /**
     * Once the permission dialog is closed
     * If all permissions are granted, we may use the app
     */
    override fun onRequestPermissionsResult(requestCode: Int, permissions: Array<out String>, grantResults: IntArray) {
        when (requestCode) {
            PERMISSION_REQUEST_CAMERA -> {
                if (grantResults.isNotEmpty() && grantResults.all { it == PackageManager.PERMISSION_GRANTED }) {
                    // Permission Granted
                    initializeCamera()
                } else {
                    // Permission Denied. Restart the app to ask again
                }
                return
            }
        }
    }

    /**
     * Calculates camera setting and updates the CameraController session
     */
    private fun calculateCameraSetting(progress: Int, key: Int) {

        // Calculate camera setting
        when (key) {
            SEEK_BAR_EXPOSURE -> {
                // camera
                if (cameraController is CameraController1) {
                    val minExp = cameraController?.cameraFeatures?.min_exposure ?: 0
                    val maxExp = cameraController?.cameraFeatures?.max_exposure ?: 0
                    exposureCompensation = minExp + progress / 100.0 * (maxExp - minExp)
                    text_view_exposure_value.text = "%.2f".format(exposureCompensation)
                }
                // camera2
                else if (cameraController is CameraController2) {
                    exposure = (100000 * exp(0.069 * progress)).toLong()
                    text_view_exposure_value.text = "%d ms".format(exposure/1000000)
                }
            }

            SEEK_BAR_FOCUS -> {
                // camera
                if (cameraController is CameraController1) {
                    // Manual focus is not supported
                    // Set focus to infinity for slider all the way to right

                    cameraController?.clearFocusAndMetering()

                    if (progress == 100) {
                        cameraController?.focusValue = "focus_mode_infinity"
                    }
                    // Otherwise set it to auto focus
                    else {
                        cameraController?.focusValue = "focus_mode_locked"
                    }

                    val params = (cameraController as CameraController1).parameters
                    val output = FloatArray(3)
                    params.getFocusDistances(output)

                    text_view_focus_value.text = "${output[Camera.Parameters.FOCUS_DISTANCE_OPTIMAL_INDEX]}"
                }
                // camera2
                else if (cameraController is CameraController2) {
                    val minFocusDist: Float = cameraController?.cameraFeatures?.minimum_focus_distance ?: 0.0f
                    // max focus distance in inf. (inputted as 0.0f)
                    // progress = 0   -> focus = minFocusDist
                    // progress = 100 -> focus = 0.0f
                    focus =  minFocusDist - ((progress / 100.0f) * minFocusDist)
                    text_view_focus_value.text = "$focus"
                }
            }

            SEEK_BAR_GAIN -> {
                // camera
                if (cameraController is CameraController1) {
                    // First we need to retrieve supported values
                    // TODO optimize?
                    val supportedGainValuesStr = cameraController?.setISO("auto")
                    val supportedGainValues: MutableList<Int> = mutableListOf()
                    var prefixPresent = false

                    // Format seems to be either ISO### or ###. So remove ISO and check if it converts to an int
                    for (gainVal: String in supportedGainValuesStr?.values!!) {
                        val formattedGainVal: String
                        val gainValInt: Int

                        if (gainVal.startsWith("ISO")) {
                            formattedGainVal = gainVal.substringAfter("ISO")
                            prefixPresent = true

                            try {
                                gainValInt = formattedGainVal.toInt()
                                Log.d(TAG, "Gain value of $gainValInt")
                                supportedGainValues.add(gainValInt)
                            } catch (e: java.lang.NumberFormatException) {
                                // Not a numbered exposure value
                                Log.d(TAG, "Non numbered gain of $formattedGainVal")
                            }
                        }
                    }

                    // If the camera is unable to change ISO, print error
                    if (supportedGainValues.size < 2) {
                        Log.e(TAG, "Error: Camera does not support changing ISO value")
                        return
                    } else {
                        val minGain = supportedGainValues.min() ?: return
                        val maxGain = supportedGainValues.max() ?: return
                        val suggestedGainValue = (minGain + (progress / 100.0) * (maxGain - minGain)).toInt()
                        var closestGainValue = supportedGainValues[0]

                        // Find closest value
                        var minDiff = Int.MAX_VALUE
                        for (gainVal: Int in supportedGainValues) {
                            if (abs(suggestedGainValue - gainVal) < minDiff) {
                                minDiff = abs(suggestedGainValue - gainVal)
                                closestGainValue = gainVal
                            }
                        }

                        gainString = if (prefixPresent) {
                            "ISO$closestGainValue"
                        } else {
                            "$closestGainValue"
                        }

                        text_view_gain_value.text = gainString
                    }
                }
                // camera2
                else if (cameraController is CameraController2) {
                    val rangeGain: Range<Int> = cameraManager.getCameraCharacteristics(cameraId)
                            .get(CameraCharacteristics.SENSOR_INFO_SENSITIVITY_RANGE)
                    val minGain = rangeGain.lower
                    val maxGain = rangeGain.upper
                    gain = (minGain + (progress / 100.0) * (maxGain - minGain)).toInt()
                    text_view_gain_value.text = "$gain"
                }
            }

            SEEK_BAR_RES -> {
                val pictureSizes = mutableListOf<CameraController.Size>()
                val previewSizes = cameraController?.cameraFeatures?.preview_sizes ?: return

                for (pictureSize in cameraController?.cameraFeatures?.picture_sizes ?: return) {
                    // Only use 16:9 resolutions (or close) that also have valid preview sizes
                    if (abs(pictureSize.width.toDouble() / pictureSize.height.toDouble() - 16.0 / 9.0) <= 0.2
                            && previewSizes.contains(pictureSize)) {
                        Log.d(TAG, "Picture size is ~16:9 $pictureSize")
                        pictureSizes.add(pictureSize)
                    }
                }

                // Choose the closest picture size
                // PictureSizes gets filled largest to smallest, so reverse here for slider
                resolution = pictureSizes[round((1.0 - progress / 100.0) * (pictureSizes.size - 1)).toInt()]

                try {
                    // camera
                    if (cameraController is CameraController1) {
                        cameraController?.stopPreview()
                        cameraController?.setPictureSize(resolution!!.width, resolution!!.height)
                        cameraController?.setPreviewSize(resolution!!.width, resolution!!.height)
                        cameraController?.startPreview()
                    }
                    // camera2
                    else if (cameraController is CameraController2) {
                        cameraController?.stopPreview()

                        if ((cameraController as CameraController2).captureSession == null) {
                            cameraController?.setPictureSize(resolution!!.width, resolution!!.height)
                            cameraController?.setPreviewSize(resolution!!.width, resolution!!.height)
                            cameraController?.startPreview()
                        }
                    }

                    text_view_res_value.text = "$resolution"
                } catch (e: CameraAccessException) {
                    e.printStackTrace()
                }
            }
        }

        // Apply the calculated setting
        try {
            // camera
            if (cameraController is CameraController1) {
                Log.d(TAG, "Using camera")
                cameraController?.setISO(gainString)
                cameraController?.exposureCompensation = exposureCompensation.toInt()
            }
            // camera2
            else if (cameraController is CameraController2) {
                Log.d(TAG, "Using camera2")
                /**
                 * CameraController will not allow auto exposure to turn off unless we also set a
                 * manual ISO first
                 */
                cameraController?.setManualISO(true, gain)
                cameraController?.exposureTime = exposure
                cameraController?.focusDistance = focus
                cameraController?.focusValue = "focus_mode_manual2"
            }

            // Calculate FOV
            val fovX = cameraController?.cameraFeatures?.view_angle_x
            val fovY = cameraController?.cameraFeatures?.view_angle_y
            val fovText = "FOV: $fovX (Horiz.) $fovY (Vert.)"
            text_view_FOV.text =  fovText
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

    /**
     * Many cameras (mostly camera but likely some camera2 hardware as well) do not have the ability
     * to record 10 second exposures. Instead, we can stack smaller exposures to reach the desired
     * exposure time (this is especially useful in situations where we do not have access to the
     * camera's current exposure time)
     */
    fun stackExposure(exposureSeconds: Double) {
        Log.d(TAG, "stackExposure()")
        val startTime = System.currentTimeMillis()
        Log.d(TAG, "bitmap")
        val width = textureView.bitmap.width
        val height = textureView.bitmap.height
        val config = textureView.bitmap.config

        // GC Issues
        val bitmapPixels = IntArray(width * height)

        // GC Issues
        //val stackedBitmapARGB = Array(width * height) { ARGBColor() }

        val stackedImage = IntArray(width * height)
        Log.d(TAG, "bitmap")


        var i = 0
        var bitmap: Bitmap

        /*while (System.currentTimeMillis() - startTime < exposureSeconds/1000.0) {*/
        while (i < 5) {
            // This makes the bitmap immutable, preventing any possible changes
            // This may not be necessary, but I'll leave it for now just in case
            bitmap = Bitmap.createBitmap(textureView.bitmap)

            // Loads pixels into stackedBitmap
            bitmap.getPixels(bitmapPixels, 0, bitmap.width, 0, 0, bitmap.width, bitmap.height)
            Log.d(TAG, "$bitmap")

            // Decode packed ARGB ints
            var color: Int
            var A: Int
            var R: Int
            var G: Int
            var B: Int

            var stackedColor: Int
            var stackedA: Int
            var stackedR: Int
            var stackedG: Int
            var stackedB: Int

            var sumA: Int
            var sumR: Int
            var sumG: Int
            var sumB: Int

            //var argbColor: ARGBColor

            for (j in bitmapPixels.indices) {
                color = bitmapPixels[j]
                stackedColor = stackedImage[j]

                // See https://developer.android.com/reference/android/graphics/Color#color-ints
                A = color shr 24 and 0xff // or color >>> 24
                R = color shr 16 and 0xff
                G = color shr 8 and 0xff
                B = color and 0xff

                stackedA = stackedColor shr 24 and 0xff // or color >>> 24
                stackedR = stackedColor shr 16 and 0xff
                stackedG = stackedColor shr 8 and 0xff
                stackedB = stackedColor and 0xff

                // Creates a new instance. DO NOT DO THIS argbColor = ARGBColor(A, R, G, B)
                // stackedBitmapARGB[j].stack(A, R, G, B)

                sumA = min(stackedA + A, 255)
                sumR = min(stackedR + R, 255)
                sumG = min(stackedG + G, 255)
                sumB = min(stackedB + B, 255)

                // Encode back into packed RGBA int
                stackedImage[j] = sumA and 0xff shl 24 or (sumR and 0xff shl 16) or (sumG and 0xff shl 8) or (sumB and 0xff)
            }

            if (cameraController is CameraController1) {

            } else if (cameraController is CameraController2) {

            }

            i++
        }

        image_view_stack.setImageBitmap(Bitmap.createBitmap(stackedImage, width, height, config))

        //cameraController?.stopPreview()
        //val canvas = textureView.lockCanvas()
        //canvas.drawARGB(255, 255, 0, 0)
    }

    //endregion

    //region old

    fun calculate(view: View) {
        val aString = editTextSemiMajor.text.toString()
        val eString = editTextEccentricity.text.toString()

        val a: Float
        val e: Float

        try {
            a = aString.toFloat()
            e = eString.toFloat()
            textViewOutput.text = computeOrbitParams(a, e)
        } catch (error: NumberFormatException) {
            val toastText = "Please enter valid numbers"
            val toastDuration = Toast.LENGTH_SHORT
            val toast = Toast.makeText(applicationContext, toastText, toastDuration)
            toast.show()
        }

    }

    /**
     * A native method that is implemented by the 'native-lib' native library,
     * which is packaged with this application.
     */
    private external fun computeOrbitParams(a: Float, e: Float): String

    //endregion

}
