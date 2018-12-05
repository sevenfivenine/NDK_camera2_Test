package com.darkfuturestudios.ndktest

import android.content.Context
import android.content.pm.PackageManager
import android.graphics.Bitmap
import android.graphics.Camera
import android.graphics.SurfaceTexture
import android.hardware.camera2.*
import android.os.*
import android.support.v4.app.ActivityCompat
import android.support.v7.app.AppCompatActivity
import android.util.Log
import android.util.Range
import android.util.Size
import android.view.Surface
import android.view.TextureView
import android.view.View
import android.widget.SeekBar
import android.widget.Toast
import com.darkfuturestudios.ndktest.CameraController.*
import kotlinx.android.synthetic.main.activity_main.*
import java.io.File
import java.io.FileOutputStream
import java.io.IOException
import java.text.SimpleDateFormat
import java.util.*
import kotlin.collections.HashMap
import kotlin.math.abs
import kotlin.math.exp
import kotlin.math.round


class MainActivity : AppCompatActivity() {

    private lateinit var stateCallback: CameraDevice.StateCallback
    private lateinit var textureView: TextureView
    private lateinit var fabTakePhoto: View
    private lateinit var seekBars: HashMap<Int, SeekBar>
    private lateinit var surfaceTextureListener: TextureView.SurfaceTextureListener
    private lateinit var cameraManager: CameraManager
    private lateinit var previewSize: Size
    private lateinit var cameraId: String

    private lateinit var galleryFolder: File

    private var cameraDevice: CameraDevice? = null
    private var cameraCaptureSession: CameraCaptureSession? = null
    private var backgroundHandler: Handler? = null
    private var backgroundThread: HandlerThread? = null
    private var captureRequestBuilder: CaptureRequest.Builder? = null
    private var captureRequest: CaptureRequest? = null
    private var cameraController: CameraController? = null
    private var hardwareSupportsCamera2: Boolean = true

    // Camera settings
    private var exposure: Long = 1000000L // camera2
    private var exposureCompensation: Double = 0.0 // camera
    private var focus: Double = 0.0
    private var gain: Int = 0 // camera2
    private var gainString: String = "" // camera
    private var resolution: CameraController.Size? = null


    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_camera)

        textureView = findViewById(R.id.texture_view)
        fabTakePhoto = findViewById(R.id.fab_take_photo)

        // Seek bars
        seekBars = hashMapOf(SEEK_BAR_EXPOSURE to findViewById(R.id.seek_bar_exposure) as SeekBar,
                SEEK_BAR_FOCUS to findViewById(R.id.seek_bar_focus) as SeekBar,
                SEEK_BAR_GAIN to findViewById(R.id.seek_bar_gain) as SeekBar,
                SEEK_BAR_RES to findViewById(R.id.seek_bar_res) as SeekBar)

        for ((key, seekBar) in seekBars) {
            seekBar.setOnSeekBarChangeListener( object : SeekBar.OnSeekBarChangeListener {
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
            /*lock()
            var outputPhoto: FileOutputStream? = null
            try {
                outputPhoto = FileOutputStream(createImageFile())
                textureView.bitmap.compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
            } catch (e: Exception) {
                e.printStackTrace()
            } finally {
                unlock()
                try {
                    outputPhoto?.close()
                } catch (e: IOException) {
                    e.printStackTrace()
                }
            }*/
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

        // This is no longer being used
        stateCallback = object : CameraDevice.StateCallback() {
            override fun onOpened(cameraDevice: CameraDevice) {
                this@MainActivity.cameraDevice = cameraDevice
                //TODO make sure we don't need this
                createPreviewSession()
        }

            override fun onDisconnected(cameraDevice: CameraDevice) {
                cameraDevice.close()
                this@MainActivity.cameraDevice = null
            }

            override fun onError(cameraDevice: CameraDevice, error: Int) {
                cameraDevice.close()
                this@MainActivity.cameraDevice = null
            }
        }
    }

    override fun onResume() {
        super.onResume()
        openBackgroundThread()
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
        //closeCameraOld()
        //closeBackgroundThread()
    }

    /**
     * Sets up camera, creates CameraController, starts preview
     */
    private fun initializeCamera() {
        val surfaceTexture = textureView.surfaceTexture

        // User has already given access
        if (checkPermissions()) {
            setUpCamera()
            createImageGallery()
            cameraController = createCameraController()
            cameraController?.setPreviewTexture(surfaceTexture)
            NDKTestUtil.setCameraDisplayOrientation(this, cameraController, hardwareSupportsCamera2)
            cameraController?.setPreviewSize(1920, 1080)
            //TODO this will crash if set to previewSize.width, previewSize.height
            cameraController?.setPictureSize(1080, 1920)
            cameraController?.startPreview()
        }
        // User has not yet given access. Prompt has been issued, waiting for callback...
        else {

        }


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
                    //TODO make this a member, give an int value
                    val cameraOpenState = "closed"
                    //TODO need this?
                    // applicationInterface.onCameraError()
                }
            }

            val useCamera2 = Build.VERSION.SDK_INT >= Build.VERSION_CODES.LOLLIPOP && hardwareSupportsCamera2

            cameraControllerLocal = if (useCamera2) {
                val previewErrorCallback = CameraController.ErrorCallback {
                    //TODO this
                    // applicationInterface.onFailedStartPreview()
                }

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
     * Use this for CameraController
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

    private fun closeCameraOld() {
        cameraCaptureSession?.close()
        cameraCaptureSession = null

        cameraDevice?.close()
        cameraDevice = null
    }

    private fun closeBackgroundThread() {
        if (backgroundHandler != null) {
            backgroundThread?.quitSafely()
            backgroundThread = null
            backgroundHandler = null
        }
    }

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
                    val streamConfigurationMap = cameraCharacteristics.get(
                            CameraCharacteristics.SCALER_STREAM_CONFIGURATION_MAP)
                    previewSize = streamConfigurationMap!!.getOutputSizes(SurfaceTexture::class.java)[0]
                    this.cameraId = cameraId
                    maxFocalLengths = focalLengths
                }

                if (!camera2Manager.allowCamera2Support(cameraId.toInt()))
                    hardwareSupportsCamera2 = false
            }

            Log.d(TAG, "Exposure time range: ${cameraManager.getCameraCharacteristics(cameraId)
                    .get(CameraCharacteristics.SENSOR_INFO_EXPOSURE_TIME_RANGE)}")
            Log.d(TAG, "Char: ${cameraManager.getCameraCharacteristics(cameraId)}")
            Log.d(TAG, "Camera: $cameraId")
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

    /**
     * Do not do this
     */
    private fun openCamera() {
        try {
            Log.d(TAG, "Checking Permissions for Camera use...")
            if (ActivityCompat.checkSelfPermission(this, android.Manifest.permission.CAMERA) == PackageManager.PERMISSION_GRANTED
            && ActivityCompat.checkSelfPermission(this, android.Manifest.permission.WRITE_EXTERNAL_STORAGE) == PackageManager.PERMISSION_GRANTED) {
                createImageGallery()
                cameraManager.openCamera(cameraId, stateCallback, backgroundHandler)
            } else {
                Log.d(TAG, "Permissions Denied!")
                ActivityCompat.requestPermissions(this@MainActivity,
                        arrayOf(android.Manifest.permission.CAMERA, android.Manifest.permission.WRITE_EXTERNAL_STORAGE), PERMISSION_REQUEST_CAMERA)
            }
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

    /**
     * Checks for permissions (camera, write external storage)
     * @return True if permissions have already been granted, false if not (will prompt user)
     */
    private fun checkPermissions(): Boolean {
        try {
            Log.d(TAG, "Checking Permissions for Camera use...")
            if (ActivityCompat.checkSelfPermission(this, android.Manifest.permission.CAMERA) == PackageManager.PERMISSION_GRANTED
                    && ActivityCompat.checkSelfPermission(this, android.Manifest.permission.WRITE_EXTERNAL_STORAGE) == PackageManager.PERMISSION_GRANTED) {
                return true
            } else {
                Log.d(TAG, "Permissions Denied!")
                ActivityCompat.requestPermissions(this@MainActivity,
                        arrayOf(android.Manifest.permission.CAMERA, android.Manifest.permission.WRITE_EXTERNAL_STORAGE), PERMISSION_REQUEST_CAMERA)
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

                // Pause preview
                if (cameraController is CameraController1) {

                } else if (cameraController is CameraController2) {
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

                // lock()
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
                    // unlock()
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

    private fun openBackgroundThread() {
        backgroundThread = HandlerThread("camera_background_thread")
        backgroundThread?.start()
        backgroundHandler = Handler(backgroundThread?.looper)
    }

    /**
     * This should be done via CameraController, so this should not be called
     */
    private fun createPreviewSession() {
        try {
            val surfaceTexture = textureView.surfaceTexture
            surfaceTexture.setDefaultBufferSize(previewSize.width, previewSize.height)
            val previewSurface = Surface(surfaceTexture)
            captureRequestBuilder = cameraDevice?.createCaptureRequest(CameraDevice.TEMPLATE_PREVIEW)
            captureRequestBuilder?.addTarget(previewSurface)

            // Adjust the exposure

            // First turn Auto Exposure off
            captureRequestBuilder?.set(CaptureRequest.CONTROL_AE_MODE, CaptureRequest.CONTROL_AE_MODE_OFF)
            Log.d(TAG, "Control mode: ${captureRequestBuilder?.get(CaptureRequest.CONTROL_AE_MODE)}")

            // Then set the exposure time (ns)
            // TODO Look into range of device exposure times
            //captureRequestBuilder?.set(CaptureRequest.SENSOR_EXPOSURE_TIME, 1*1000000000)
            captureRequestBuilder?.set(CaptureRequest.SENSOR_EXPOSURE_TIME, exposure)
            //captureRequestBuilder?.set(CaptureRequest.SENSOR_FRAME_DURATION, 100000000)

            cameraDevice?.createCaptureSession(Collections.singletonList(previewSurface),
                    object : CameraCaptureSession.StateCallback() {

                        override fun onConfigured(cameraCaptureSession: CameraCaptureSession) {
                            if (cameraDevice ==
                                    null) {
                                return
                            }

                            try {
                                captureRequest = captureRequestBuilder?.build()
                                this@MainActivity.cameraCaptureSession = cameraCaptureSession
                                this@MainActivity.cameraCaptureSession?.setRepeatingRequest(captureRequest, null, backgroundHandler)
                            } catch (e: CameraAccessException) {
                                e.printStackTrace()
                            }

                        }

                        override fun onConfigureFailed(cameraCaptureSession: CameraCaptureSession) {

                        }
                    }, backgroundHandler)
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

    override fun onRequestPermissionsResult(requestCode: Int, permissions: Array<out String>, grantResults: IntArray) {
        when (requestCode) {
            PERMISSION_REQUEST_CAMERA -> {
                if (grantResults.isNotEmpty() && grantResults.all { it == PackageManager.PERMISSION_GRANTED } ) {
                    // Permission Granted
                    initializeCamera()
                    //if(textureView.isAvailable) openCamera()
                } else {
                    // Permission Denied
                }
                return
            }
        }
    }

    /**
     * Creates a gallery in the phones Pictures directory for this app
     */
    private fun createImageGallery() {
        // On any Android phone, there is a public directory for pictures
        val storageDirectory = Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_PICTURES)
        galleryFolder = File(storageDirectory, resources.getString(R.string.app_name))
        if (!galleryFolder.exists()) {
            val wasCreated = galleryFolder.mkdirs()
            if(!wasCreated) {
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
     * When a picture is taken, the TextureView is frozen for a brief time to indicate that the
     * picture has been captured
     */
    private fun lock() {
        try {
            // null 2nd argument indicates we do not need to use image metadata
            cameraCaptureSession?.capture(captureRequestBuilder?.build(), null, backgroundHandler)
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }
    }

    /**
     * After a brief pause, the TextureView unfreezes and repeatedly updates, indicating the user
     * is in capture mode again
     */
    private fun unlock() {
        try {
            // null 2nd argument indicates we do not need to use image metadata
            cameraCaptureSession?.setRepeatingRequest(captureRequestBuilder?.build(), null, backgroundHandler)
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }
    }

    /**
     * Calculates camera setting and updates the current session
     */
    private fun calculateCameraSetting(progress: Int, key: Int) {
        Log.d(TAG, "Progress: " + progress)

        // Calculate camera setting
        when (key) {
            SEEK_BAR_EXPOSURE -> {
                Log.d(TAG, "Exposure: $exposure")

                // camera
                if (cameraController is CameraController1) {
                    //TODO what if these are null?
                    val minExp = cameraController?.cameraFeatures?.min_exposure ?: 0
                    val maxExp = cameraController?.cameraFeatures?.max_exposure ?: 0
                    exposureCompensation = minExp + progress/100.0*(maxExp - minExp)
                    Log.d(TAG, "Min ExpComp $minExp Max ExpComp $maxExp")
                }
                // camera2
                else if (cameraController is CameraController2) {
                    exposure = (100000* exp(0.069*progress)).toLong()
                }
            }
            SEEK_BAR_FOCUS -> {
                //focus = progress
            }
            SEEK_BAR_GAIN -> {
                // camera
                if (cameraController is CameraController1) {
                    // First we need to retrieve supported values
                    // TODO optimize?
                    val supportedGainValuesStr = cameraController?.setISO("auto")
                    val supportedGainValues: MutableList<Int> = mutableListOf()
                    var prefixPresent = false

                    // Format seems to be either ISO400 or 400. So remove ISO and check if it converts to an int
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
                        val suggestedGainValue = (minGain + (progress/100.0)*(maxGain - minGain)).toInt()
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
                    }
                }
                // camera2
                else if (cameraController is CameraController2) {
                    val rangeGain: Range<Int> = cameraManager.getCameraCharacteristics(cameraId).get(CameraCharacteristics.SENSOR_INFO_SENSITIVITY_RANGE)
                    val minGain = rangeGain.lower
                    val maxGain = rangeGain.upper
                    gain = (minGain + (progress/100.0)*(maxGain - minGain)).toInt()
                }
            }
            SEEK_BAR_RES -> {
                val pictureSizes = mutableListOf<CameraController.Size>()
                val previewSizes = cameraController?.cameraFeatures?.preview_sizes ?: return

                for (pictureSize in cameraController?.cameraFeatures?.picture_sizes ?: return) {
                    // Only use 16:9 resolutions (or close) that also have valid preview sizes
                    if (abs(pictureSize.width.toDouble()/pictureSize.height.toDouble() - 16.0/9.0) <= 0.2
                    && previewSizes.contains(pictureSize)) {
                        Log.d(TAG, "Picture size is ~16:9 $pictureSize")
                        pictureSizes.add(pictureSize)
                    }
                }

                // Choose the closest picture size
                // PictureSizes gets filled largest to smallest, so reverse here for slider
                resolution = pictureSizes[round((1.0 - progress/100.0)*(pictureSizes.size - 1)).toInt()]

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
                        }

                        if ((cameraController as CameraController2).captureSession == null) {
                            cameraController?.startPreview()
                        }
                    }

                } catch (e: CameraAccessException) {
                    e.printStackTrace()
                }
            }
        }

        // Apply the calculated setting
        try {
            //val surfaceTexture = textureView.surfaceTexture ?: return
            //val previewSurface = Surface(surfaceTexture)
            //captureRequestBuilder = cameraDevice?.createCaptureRequest(CameraDevice.TEMPLATE_PREVIEW)
            //captureRequestBuilder?.addTarget(previewSurface)

            //captureRequestBuilder?.set(CaptureRequest.CONTROL_AE_MODE, CaptureRequest.CONTROL_AE_MODE_OFF)
            //captureRequestBuilder?.set(CaptureRequest.CONTROL_MODE, CaptureRequest.CONTROL_MODE_OFF)
            //captureRequestBuilder?.set(CaptureRequest.SENSOR_EXPOSURE_TIME, exposure)

            // camera
            if (cameraController is CameraController1) {
                Log.d(TAG, "Using camera")
                val isoVal = "ISO1600"
                //cameraController?.autoExposureLock = false
                cameraController?.setISO(gainString)
                cameraController?.exposureCompensation = exposureCompensation.toInt()
                //cameraController?.autoExposureLock = true
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
            }

            Log.d(TAG, "EXP TIME: ${captureRequestBuilder?.get(CaptureRequest.SENSOR_EXPOSURE_TIME)}")
            Log.d(TAG, "Control mode: ${captureRequestBuilder?.get(CaptureRequest.CONTROL_AE_MODE)}")

            //captureRequest = captureRequestBuilder?.build()
            //Log.d(TAG, "Request: $captureRequest")
            //cameraCaptureSession?.setRepeatingRequest(captureRequest, null, backgroundHandler)
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

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

    companion object {

        const val TAG: String = "MainActivity"
        const val PERMISSION_REQUEST_CAMERA: Int = 7000

        const val SEEK_BAR_EXPOSURE = 0
        const val SEEK_BAR_FOCUS = 1
        const val SEEK_BAR_GAIN = 2
        const val SEEK_BAR_RES = 3

        // Used to load the 'native-lib' library on application startup.
        init {
            System.loadLibrary("native-lib")
        }

    }

}
