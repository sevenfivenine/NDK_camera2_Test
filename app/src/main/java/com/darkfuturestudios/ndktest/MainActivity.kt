package com.darkfuturestudios.ndktest

import android.content.Context
import android.content.pm.PackageManager
import android.content.res.ColorStateList
import android.graphics.Bitmap
import android.graphics.SurfaceTexture
import android.hardware.Camera
import android.hardware.camera2.*
import android.os.*
import android.support.design.widget.FloatingActionButton
import android.support.v4.app.ActivityCompat
import android.support.v7.app.AppCompatActivity
import android.util.Log
import android.util.Range
import android.view.Menu
import android.view.MenuItem
import android.view.TextureView
import android.view.View
import android.widget.SeekBar
import android.widget.Switch
import android.widget.Toast
import com.darkfuturestudios.ndktest.CameraController.*
import kotlinx.android.synthetic.main.activity_camera.*
import kotlinx.android.synthetic.main.activity_main.*
import java.io.File
import java.io.FileOutputStream
import java.io.IOException
import java.text.SimpleDateFormat
import java.util.*
import kotlin.collections.HashMap
import kotlin.math.*

class MainActivity : AppCompatActivity() {

    //region companion

    companion object {

        const val TAG: String = "MainActivity"
        const val PERMISSION_REQUEST_CAMERA: Int = 7000

        const val SEEK_BAR_EXPOSURE = 0
        const val SEEK_BAR_FOCUS = 1
        const val SEEK_BAR_GAIN = 2
        const val SEEK_BAR_RES = 3

        // For experimenting
        // 0: Nothing
        // 1: Make the camera oscillate back and forth so the motion of the motor is visible
        const val EASTER_EGG = 0

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
    private lateinit var autoSwitches: HashMap<Int, Switch>
    private lateinit var surfaceTextureListener: TextureView.SurfaceTextureListener
    private lateinit var cameraManager: CameraManager
    private lateinit var cameraId: String
    private lateinit var galleryFolder: File

    private var cameraController: CameraController? = null
    private var hardwareSupportsCamera2: Boolean = true

    // Camera settings
    private var settingsVisible: Boolean = true // Toggles the UI (sliders, switches, labels)
    private var exposure: Long = 1000000L // camera2; Exposure time (Nanoseconds)
    private var exposureCompensation: Double = 0.0 // camera
    private var focus: Float = 0.0f
    private var gain: Int = 0 // camera2
    private var gainString: String = "" // camera
    private var resolution: CameraController.Size? = null
    private var autoModes: HashMap<Int, Boolean>? = null
    private var fovX: Float? = 0.0f  // camera angular fieid of view width, in degrees
    private var fovY: Float? = 0.0f  // camera angular field of view height, in degrees

    // Stacking

    /**
     * How long to stack exposure (Nanoseconds)
     */
    private var stackDuration: Long = 0L

    /**
     * Time when stacking started (Milliseconds)
     */
    private var stackingStartTime: Long = 0L

    /**
     * Number of individual exposures required to generate total stacked exposure
     */
    private var exposuresNeeded: Int = 0

    /**
     * Number of individual exposures actually taken so far in the current stacking sequence
     */
    private var exposuresTaken: Int = 0

    /**
     * Amount of time the shutter stays open for each stack capture (Nanoseconds)
     */
    private var stackingExposureTime: Long = 0L

    /**
     * When stacking is active, the image data is stored here
     */
    private var stackedImage: IntArray = IntArray(0)

    /**
     * When true, FAB becomes stack button
     * When false, FAB becomes photo button
     */
    private var stackMode: Boolean = false

    /**
     * Schedules stack frames
     */
    @Deprecated(message = "Use handler") private var timer: Timer? = null

    /**
     * Use this for scheduling instead
     */
    private var handler: Handler? = null

    /**
     * Units: milliseconds
     */
    private var lastStackFrameTimestamp: Long = 0L

    //endregion

    //region lifecycle

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_camera)
        setSupportActionBar(findViewById(R.id.toolbar_camera))

        // Find views
        textureView = findViewById(R.id.texture_view)
        fabTakePhoto = findViewById(R.id.fab_take_photo)

        // Create seekBars hash map
        seekBars = hashMapOf(SEEK_BAR_EXPOSURE to seek_bar_exposure,
                SEEK_BAR_FOCUS to seek_bar_focus,
                SEEK_BAR_GAIN to seek_bar_gain,
                SEEK_BAR_RES to seek_bar_res)

        // Create auto switch hash map
        autoSwitches = hashMapOf(SEEK_BAR_EXPOSURE to switch_auto_exposure,
                SEEK_BAR_FOCUS to switch_auto_focus,
                SEEK_BAR_GAIN to switch_auto_gain,
                SEEK_BAR_RES to switch_auto_resolution)

        // Create auto modes hash map
        autoModes = hashMapOf(SEEK_BAR_EXPOSURE to false,
                SEEK_BAR_FOCUS to false,
                SEEK_BAR_GAIN to false,
                SEEK_BAR_RES to false)

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

        for ((key, autoSwitch) in autoSwitches) {
            // Set listener
            autoSwitch.setOnCheckedChangeListener { compoundButton, isChecked ->
                Log.d(TAG, "Switch $key is now $isChecked")
                if (isChecked) {
                    autoModes!![key] = true

                    // Auto exposure and gain must have the same state
                    if (key == SEEK_BAR_EXPOSURE) autoModes!![SEEK_BAR_GAIN] = true
                    else if (key == SEEK_BAR_GAIN) autoModes!![SEEK_BAR_EXPOSURE] = true

                    calculateCameraSetting(seekBars[key]!!.progress, key)

                    // Disable manual
                    seekBars[key]!!.isEnabled = false
                } else {
                    autoModes!![key] = false

                    // Auto exposure and gain must have the same state
                    if (key == SEEK_BAR_EXPOSURE) autoModes!![SEEK_BAR_GAIN] = false
                    else if (key == SEEK_BAR_GAIN) autoModes!![SEEK_BAR_EXPOSURE] = false

                    calculateCameraSetting(seekBars[key]!!.progress, key)

                    // Enable manual
                    seekBars[key]!!.isEnabled = true
                }
            }
        }

        handler = Handler()

        fabTakePhoto.setOnClickListener {
            // First make the button turn red

            if (stackMode) {
                stackExposure()
            } else {
                fabTakePhoto.backgroundTintList = ColorStateList.valueOf(resources.getColor(R.color.colorStacking))

                handler?.postDelayed({
                    takePhoto(true, false)
                }, 1000)
            }

            // Wait one second for the device to stabilize after user touches it


        }

        /*fab_stack.setOnClickListener {
            stackExposure()
        }*/

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

    override fun onCreateOptionsMenu(menu: Menu?): Boolean {
        menuInflater.inflate(R.menu.menu_camera, menu)
        return super.onCreateOptionsMenu(menu)
    }

    override fun onOptionsItemSelected(item: MenuItem?): Boolean {
        if (item?.itemId == R.id.action_settings_visibility) {
            settingsVisible = !settingsVisible

            toggleUI(if (settingsVisible) View.VISIBLE else View.INVISIBLE)
        }

        return super.onOptionsItemSelected(item)
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
     * Toggles the UI (sliders, switches, labels for changing camera settings)
     */
    private fun toggleUI(visibility: Int) {
        for ((key, seekBar) in seekBars) {
            seekBar.visibility = visibility
        }

        for ((key, autoSwitch) in autoSwitches) {
            autoSwitch.visibility = visibility
        }

        text_view_exposure.visibility = visibility
        text_view_exposure_value.visibility = visibility
        text_view_focus.visibility = visibility
        text_view_focus_value.visibility = visibility
        text_view_gain.visibility = visibility
        text_view_gain_value.visibility = visibility
        text_view_res.visibility = visibility
        text_view_res_value.visibility = visibility
        text_view_FOV.visibility = visibility
        switch_hide_preview.visibility = visibility
    }

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

            // Add indicator to Action Bar that shows which camera API is used for the device
            toolbar_camera.title = if (cameraController is CameraController1) "camera (old) API"
            else if (cameraController is CameraController2) "camera2 API" else "Error: no camera?"

            // Oscillate camera
            if (EASTER_EGG == 1 && cameraController is CameraController2) {
                Log.d(TAG, "Oscillate!")
                Thread(Runnable {
                    val timer = Timer()

                    timer.scheduleAtFixedRate(object : TimerTask() {
                        override fun run() {
                            NDKTestUtil.LOG = false
                            val minFocusDist: Float = cameraController?.cameraFeatures?.minimum_focus_distance ?: 0.0f
                            focus = ((minFocusDist/2) + (minFocusDist/2) * sin(6 * System.currentTimeMillis()/1000.0)).toFloat()
                            seekBars[SEEK_BAR_FOCUS]?.progress = (100.0f * (minFocusDist - focus)/minFocusDist).toInt()
                            //text_view_focus_value.text = "$focus"
                            cameraController?.focusDistance = focus
                            cameraController?.focusValue = "focus_mode_manual2"
                            Log.d(TAG, "Oscillating! $focus")
                            NDKTestUtil.LOG = true
                        }

                    }, 0, 100)
                }).start()
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
    private fun takePhoto(save: Boolean, stack: Boolean) {
        Log.d(TAG, "Take photo")

        if (!textureView.isAvailable) {
            Log.d(TAG, "Texture view not yet available")
            return
        }

        // Only needed in photo mode. When stacking simply read data directly from TextureView
        if (!stack) {
            val pictureCallback = object : CameraController.PictureCallback {
                override fun onStarted() {
                    Log.d(TAG, "PictureCallback.onStarted()")
                }

                override fun onCompleted() {
                    Log.d(TAG, "PictureCallback.onCompleted()")

                    // For stacking. When finished taking photo, do next capture
                    /*if (stack) {
                        runSingleStack()
                    }*/
/*
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
*/
                }

                /**
                 * @param data contains EXIF data from the camera, including exposure time
                 */
                override fun onPictureTaken(data: ByteArray?) {
                    Log.d(TAG, "PictureCallback.onPictureTaken()")

                    if (save) {
                        var outputPhoto: FileOutputStream? = null
                        try {
                            outputPhoto = FileOutputStream(createImageFile())
                            /** OpenCamera uses a much more sophisticated method of saving images
                             *  This is much simpler, but less versatile
                             *  We don't even use the data argument, just capture from textureView instead
                             */
                            //textureView.bitmap.compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
                            textureView.getBitmap(resolution!!.height, resolution!!.width).compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
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
                }

                override fun onRawPictureTaken(raw_image: RawImage?) {}

                override fun onBurstPictureTaken(images: MutableList<ByteArray>?) {}

                override fun onFrontScreenTurnOn() {}

            }

            val errorCallback = CameraController.ErrorCallback { Log.e(TAG, "Error from takePicture()") }

            cameraController?.takePicture(pictureCallback, errorCallback)
        }
        // Get data from TextureView
        else {
            if (save) {
                var outputPhoto: FileOutputStream? = null
                try {
                    outputPhoto = FileOutputStream(createImageFile())
                    /** OpenCamera uses a much more sophisticated method of saving images
                     *  This is much simpler, but less versatile
                     *  We don't even use the data argument, just capture from textureView instead
                     */
                    //textureView.bitmap.compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
                    textureView.getBitmap(resolution!!.height, resolution!!.width).compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
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

            /**
             * For stacking: Data has been captured. Now send it to C
             */
            if (stack) {
                Log.d(TAG, "Processing single stack frame")
                val previewBitmap = textureView.getBitmap(resolution!!.height, resolution!!.width)
                val width = previewBitmap.width
                val height = previewBitmap.height
                val config = previewBitmap.config

                // This makes the bitmap immutable, preventing any possible changes
                // This may not be necessary, but I'll leave it for now just in case
                val bitmap = Bitmap.createBitmap(textureView.getBitmap(resolution!!.height, resolution!!.width))

                val bitmapPixels = IntArray(width * height)

                // Loads pixels into stackedBitmap
                bitmap.getPixels(bitmapPixels, 0, bitmap.width, 0, 0, bitmap.width, bitmap.height)
                Log.d(TAG, "$bitmap")

                stackImageBuffers( bitmapPixels, width, height, stackedImage )

                // Check if stacking is completed
                if (exposuresTaken >= exposuresNeeded) {
                    Log.d(TAG,"Stacking complete")
                    processStackedImage(width, height, config)
                }
            }
        }

        if (!stack)
        // Change FAB color back to blue
        fabTakePhoto.backgroundTintList = ColorStateList.valueOf(resources.getColor(R.color.colorPrimary))
    }

    private fun processStackedImage(width: Int, height: Int, config: Bitmap.Config) {
        Log.d(TAG, "processStackImage()")
        // Process stacked image in native code
        var imageFile: File = createImageFile()
        val imagePath: String = imageFile.absolutePath
        val galleryPath: String = galleryFolder.absolutePath
        processStackedImage ( stackedImage, width, height, fovX ?: 0.0f, fovY ?: 0.0f, stackDuration/1.0e9f, gain, imagePath, galleryPath )

        // Display stacked image

        val stackedBitmap = Bitmap.createBitmap(stackedImage, width, height, config)
        image_view_stack.setImageBitmap(stackedBitmap)
        switch_hide_preview.isChecked = true

        // Save stacked image

        var outputPhoto: FileOutputStream? = null
        try {
            outputPhoto = FileOutputStream(imageFile)
            /** OpenCamera uses a much more sophisticated method of saving images
             *  This is much simpler, but less versatile
             *  We don't even use the data argument, just capture from textureView instead
             */
            stackedBitmap.compress(Bitmap.CompressFormat.JPEG, 100, outputPhoto)
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
                // Auto exposure
                if (autoModes!![SEEK_BAR_EXPOSURE] == true) {
                    text_view_exposure_value.text = "Auto"
                }

                // Manual exposure
                else {
                    val minExposureTime = cameraController?.cameraFeatures?.min_exposure_time ?: 0
                    val maxExposureTime = cameraController?.cameraFeatures?.max_exposure_time ?: 0

                    // camera
                    if (cameraController is CameraController1) {
                        val minExp = cameraController?.cameraFeatures?.min_exposure ?: 0
                        val maxExp = cameraController?.cameraFeatures?.max_exposure ?: 0
                        exposureCompensation = minExp + progress / 100.0 * (maxExp - minExp)
                        text_view_exposure_value.text = "%.2f".format(exposureCompensation)
                    }
                    // camera2
                    else if (cameraController is CameraController2) {
                        /*

                        Not true! See below!

                        // If max exposure time of the camera is greater than 10 seconds
                        // We do not ever need to use stacking
                        if (maxExposureTime >= 10000000000) {
                            // linear: exposure = (minExposureTime + (progress / 100.0f) * (10000000000 - minExposureTime)).toLong()
                            // exponential:
                            exposure = (minExposureTime * exp(0.01*ln(10000000000.0/minExposureTime) * progress)).toLong()
                        }

                        // If max exposure time of the camera is less than 10 seconds
                        // Stacking is necessary past max exposure time
                        else {
                            exposure = (100000 * exp(0.069 * progress)).toLong()
                        }*/

                        // We always need stacking because, at least on some devices, showing the
                        // preview at 10 seconds will cause a camera firmware crash
                        // Since we want to support a variety of Android devices, it is in our best
                        // interest to go the safe route and stack

                        // The point where the slider switches from changing exposure time to stacking
                        val minStackingProgress = 50.1

                        // 1/2 a second, or the max exposure time of the device
                        val previewMaxExposureTime = min((1.0/2) * 1000000000L, maxExposureTime.toDouble())

                        // Use 1/2 a second, or the max exposure time of the device
                        stackingExposureTime = previewMaxExposureTime.toLong()

                        val maxStackTime = 4L * 1000000000L // 4 seconds

                        /**
                         * Change exposure time
                         * This portion of the slider goes from minExposureTime to
                         * previewMaxExposureTime (1/12 second is a good value, used on OpenCamera)
                         */
                        if (progress < minStackingProgress) {
                            stackMode = false
                            (fabTakePhoto as FloatingActionButton).setImageResource(R.drawable.ic_camera_light)
                            // Not stacking
                            stackDuration = 0L
                            // Progress through the exposure time portion of the seek bar
                            // 0 to 1
                            val exposureTimeProgress = progress/minStackingProgress
                            exposure = (minExposureTime * exp(ln(previewMaxExposureTime/minExposureTime)
                                    * exposureTimeProgress)).toLong()

                            text_view_exposure_value.text = "%d ms".format(exposure/1000000)
                        }

                        /**
                         * Stack exposure
                         */
                        else {
                            stackMode = true
                            (fabTakePhoto as FloatingActionButton).setImageResource(R.drawable.ic_stack_light)
                            // Set the exposure time (for each capture)
                            exposure = stackingExposureTime
                            // Progress through the stacking portion of the seek bar
                            // 0 to 1
                            val stackProgress = (progress - minStackingProgress)/(100.0f - minStackingProgress)

                            // Linear
                            stackDuration = (previewMaxExposureTime + (stackProgress * (maxStackTime - previewMaxExposureTime))).toLong()

                            text_view_exposure_value.text = "%.2f ms".format(stackDuration/1000000.0)
                        }
                    }
                }
            }

            SEEK_BAR_FOCUS -> {
                // Auto focus
                if (autoModes!![SEEK_BAR_FOCUS] == true) {
                    text_view_focus_value.text = "Auto"
                }

                // Manual focus
                else {
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
                        text_view_focus_value.text = "%.2f m".format(1.0/focus)
                    }
                }
            }

            SEEK_BAR_GAIN -> {
                // Auto gain
                if (autoModes!![SEEK_BAR_GAIN] == true) {
                    text_view_gain_value.text = "Auto"
                }

                // Manual gain
                else {
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
                            } else {
                                formattedGainVal = gainVal
                                prefixPresent = false
                            }

                            try {
                                gainValInt = formattedGainVal.toInt()
                                Log.d(TAG, "Gain value of $gainValInt")
                                supportedGainValues.add(gainValInt)
                            } catch (e: java.lang.NumberFormatException) {
                                // Not a numbered exposure value
                                Log.d(TAG, "Non numbered gain of $formattedGainVal")
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
            }

            SEEK_BAR_RES -> {
                // Auto resolution
                if (autoModes!![SEEK_BAR_RES] == true) {
                    // Use largest possible resolution (0th element in array)
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

                    resolution = pictureSizes[0]

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

                // Manual resolution
                else {
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
        }

        // Apply the calculated setting
        try {
            // camera
            if (cameraController is CameraController1) {
                Log.d(TAG, "Using camera")

                // Auto exposure/gain (just set exposure compensation to 0 EV)
                if (autoModes!![SEEK_BAR_EXPOSURE] == true || autoModes!![SEEK_BAR_GAIN] == true) {
                    cameraController?.exposureCompensation = 0
                    cameraController?.setISO("auto")

                    // These are coupled, so make sure both are checked
                    autoSwitches[SEEK_BAR_EXPOSURE]!!.isChecked = true
                    autoSwitches[SEEK_BAR_GAIN]!!.isChecked = true
                }

                // Manual exposure
                else {
                    cameraController?.exposureCompensation = exposureCompensation.toInt()
                    cameraController?.setISO(gainString)

                    // These are coupled, so make sure they are both unchecked
                    autoSwitches[SEEK_BAR_EXPOSURE]!!.isChecked = false
                    autoSwitches[SEEK_BAR_GAIN]!!.isChecked = false
                }

                // Auto focus
                if (autoModes!![SEEK_BAR_FOCUS] == true) {

                }

                // Manual focus
                else {

                }
            }

            // camera2
            else if (cameraController is CameraController2) {
                Log.d(TAG, "Using camera2")
                Log.d(TAG, "Min exposure time: ${cameraController?.cameraFeatures?.min_exposure_time}")
                Log.d(TAG, "Max exposure time: ${cameraController?.cameraFeatures?.max_exposure_time}")
                Log.d(TAG, "Current exposure time: ${cameraController?.exposureTime}")

                // Auto exposure/gain
                if (autoModes!![SEEK_BAR_EXPOSURE] == true || autoModes!![SEEK_BAR_GAIN] == true) {
                    // Use this to turn on auto exposure and auto ISO
                    cameraController?.setManualISO(false, 0)

                    // These are coupled, so make sure both are checked
                    autoSwitches[SEEK_BAR_EXPOSURE]!!.isChecked = true
                    autoSwitches[SEEK_BAR_GAIN]!!.isChecked = true
                }

                // Manual exposure/gain
                else {
                    /**
                     * CameraController will not allow auto exposure to turn off unless we also set a
                     * manual ISO first
                     */
                    cameraController?.setManualISO(true, gain)
                    cameraController?.exposureTime = exposure

                    // These are coupled, so make sure they are both unchecked
                    autoSwitches[SEEK_BAR_EXPOSURE]!!.isChecked = false
                    autoSwitches[SEEK_BAR_GAIN]!!.isChecked = false
                }

                // Auto focus
                if (autoModes!![SEEK_BAR_FOCUS] == true) {
                    cameraController?.focusValue = "focus_mode_auto"
                }

                // Manual focus
                else {
                    cameraController?.focusDistance = focus
                    cameraController?.focusValue = "focus_mode_manual2"
                }
            }

            // Calculate FOV
            fovX = cameraController?.cameraFeatures?.view_angle_y
            fovY = cameraController?.cameraFeatures?.view_angle_x
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
    private fun stackExposure() {
        Log.d(TAG, "stackExposure()")
        //val startTime = System.currentTimeMillis()

        val previewBitmap = textureView.getBitmap(resolution!!.height, resolution!!.width)
        val width = previewBitmap.width
        val height = previewBitmap.height
        val config = previewBitmap.config

        //Log.d(TAG, "Resolution: width ${resolution?.width} height ${resolution?.height}, bitmap res width $width x height $height")
        //Log.d(TAG, "FOV width $fovX height $fovY")

        // GC Issues
        //val bitmapPixels = IntArray(width * height)

        // Initialize stacked image member to store data
        stackedImage = IntArray(width * height)

        //var bitmap: Bitmap

        stackingStartTime = System.currentTimeMillis()
        exposuresNeeded = ( ( stackDuration + stackingExposureTime - 1 ) / stackingExposureTime ).toInt()  // round up
        exposuresTaken = 0

        // Will cancel when done (see runSingleStack() )
        // Max of 100 attempts
        for (i in 0..100) {
            handler?.postDelayed({
                Log.d(TAG, "Stack frame scheduled")
                runSingleStack()
            }, (stackingExposureTime*i)/1000000)
        }



        /*timer = Timer()
        timer?.scheduleAtFixedRate(object : TimerTask() {
            override fun run() {
                Log.d(TAG, "Not on UI thread")
                runOnUiThread {
                    Log.d(TAG, "Stack frame scheduled")
                    runSingleStack()
                }
            }
        }, 0, stackingExposureTime/1000000)*/

        //runSingleStack()

        //cameraController?.stopPreview()
        //val canvas = textureView.lockCanvas()
        //canvas.drawARGB(255, 255, 0, 0)
    }

    /**
     * Runs one iteration of stacking
     * Or, if stacking is finished: Processes, displays, and saves stacked image
     * This calls takePhoto()
     * First, we must check if the stacking is finished
     */
    fun runSingleStack() {
        Log.d(TAG, "Exposure time: ${cameraController?.exposureTime}")
        // Still stacking?
        if (exposuresTaken < exposuresNeeded) {
            val acceptableError = 100 // ms
            val delay = abs((System.currentTimeMillis() - lastStackFrameTimestamp) - stackingExposureTime/1000000)

            // Once stacking delay has "stabilized", we can actually start
            if (delay <= acceptableError) {
                fabTakePhoto.backgroundTintList = ColorStateList.valueOf(resources.getColor(R.color.colorStacking))

                // For final frame in stack, set exposure time as needed to complete the whole stack to its intended duration.
                if ( exposuresTaken == exposuresNeeded - 1 )
                    cameraController?.exposureTime = stackDuration - exposuresTaken * stackingExposureTime

                Log.d(TAG, "Acceptable Stack delay $delay")
                Log.d(TAG, "Time diff ${System.currentTimeMillis() - lastStackFrameTimestamp}")
                Log.d(TAG, "Stacking exposure " + ( exposuresTaken + 1 ) + " of " + exposuresNeeded + ": " + cameraController?.exposureTime + " nanosec" )

                exposuresTaken++
                lastStackFrameTimestamp = System.currentTimeMillis()

                takePhoto(false, true)
            }

            // But if the delay is still to big, don't actually stack yet
            else {
                Log.d(TAG, "Too Much Stack delay! $delay")
                Log.d(TAG, "Current time ${System.currentTimeMillis()}")
                Log.d(TAG, "Last stack frame timestamp $lastStackFrameTimestamp")
                Log.d(TAG, "Time diff ${System.currentTimeMillis() - lastStackFrameTimestamp}")

                // Set the button to waiting color
                fabTakePhoto.backgroundTintList = ColorStateList.valueOf(resources.getColor(R.color.colorWaiting))
                lastStackFrameTimestamp = System.currentTimeMillis()
            }


        }

        // Stacking completed
        // BUT the stacked image is probably not done being captured!
        else {

            // Set camera exposure time back to its expected value

            cameraController?.exposureTime = exposure

            // NOT necessarily ready to process yet!!!

            /*val previewBitmap = textureView.getBitmap(resolution!!.height, resolution!!.width)
            val width = previewBitmap.width
            val height = previewBitmap.height
            val config = previewBitmap.config

            var imageFile: File = createImageFile()
            val imagePath: String = imageFile.absolutePath
            val galleryPath: String = galleryFolder.absolutePath*/
            //processStackedImage ( stackedImage, width, height, fovX ?: 0.0f, fovY ?: 0.0f, stackDuration/1.0e9f, gain, imagePath, galleryPath )

            // Stop the timer to stop stacking
            handler?.removeCallbacksAndMessages(null)

            // Change FAB color back to blue
            fabTakePhoto.backgroundTintList = ColorStateList.valueOf(resources.getColor(R.color.colorPrimary))
        }

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

    private external fun processStackedImage (data:IntArray, width:Int, height:Int, widthAngle:Float, heightAngle:Float, exposureSeconds:Float, iso:Int, imageFilePath:String, logDir:String ): Int

    private external fun stackImageBuffers (data:IntArray, width:Int, height:Int, stackedData:IntArray ): Int

    //endregion
}
