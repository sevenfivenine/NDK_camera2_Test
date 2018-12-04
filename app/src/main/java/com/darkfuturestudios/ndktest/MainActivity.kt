package com.darkfuturestudios.ndktest

import android.content.Context
import android.content.pm.PackageManager
import android.graphics.Bitmap
import android.graphics.SurfaceTexture
import android.hardware.camera2.*
import android.os.*
import android.support.v4.app.ActivityCompat
import android.support.v7.app.AppCompatActivity
import android.util.Log
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
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.log


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
    private var exposure: Long = 1000000L
    private var focus: Double = 0.0
    private var gain: Double = 0.0
    private var res: Double = 0.0


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
            lock()
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
            }
        }

        cameraManager = applicationContext.getSystemService(Context.CAMERA_SERVICE) as CameraManager

        surfaceTextureListener = object : TextureView.SurfaceTextureListener {
            override fun onSurfaceTextureAvailable(surfaceTexture: SurfaceTexture, width: Int, height: Int) {
                setUpCamera()
                cameraController = createCameraController()
                cameraController?.setPreviewTexture(surfaceTexture)
                cameraController?.setPreviewSize(previewSize.width, previewSize.height)
                //TODO this will crash if set to previewSize.width, previewSize.height
                cameraController?.setPictureSize(1080, 1920)
                cameraController?.startPreview()
                //openCamera()
            }

            override fun onSurfaceTextureSizeChanged(surfaceTexture: SurfaceTexture, width: Int, height: Int) {

            }

            override fun onSurfaceTextureDestroyed(surfaceTexture: SurfaceTexture): Boolean {
                return false
            }

            override fun onSurfaceTextureUpdated(surfaceTexture: SurfaceTexture) {

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
            val surfaceTexture = textureView.surfaceTexture

            setUpCamera()
            cameraController = createCameraController()
            cameraController?.setPreviewTexture(surfaceTexture)
            cameraController?.setPreviewSize(previewSize.width, previewSize.height)
            //TODO this will crash if set to previewSize.width, previewSize.height
            cameraController?.setPictureSize(1080, 1920)
            cameraController?.startPreview()
            //openCamera()
        } else {
            textureView.surfaceTextureListener = surfaceTextureListener
        }
    }

    override fun onStop() {
        super.onStop()
        closeCamera()
        closeBackgroundThread()
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

                CameraController2(applicationContext, cameraId.toInt(), previewErrorCallback, cameraErrorCallback)
            } else {
                CameraController1(cameraId.toInt(), cameraErrorCallback)
            }
        } catch (e: CameraControllerException) {
            e.printStackTrace()
            cameraControllerLocal = null
        }

        return cameraControllerLocal
    }

    private fun closeCamera() {
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
     * TODO move permission checks elsewhere
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
                    if(textureView.isAvailable) openCamera()
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

        when (key) {
            SEEK_BAR_EXPOSURE -> {
                //
                exposure = (100000* exp(0.069*progress)).toLong()
                Log.d(TAG, "Exposure: $exposure")
            }
            SEEK_BAR_FOCUS -> {
                //focus = progress
            }
            SEEK_BAR_GAIN -> {
                //gain =
            }
            SEEK_BAR_RES -> {
                //res = progress
            }
        }

        try {
            val surfaceTexture = textureView.surfaceTexture ?: return
            val previewSurface = Surface(surfaceTexture)
            captureRequestBuilder = cameraDevice?.createCaptureRequest(CameraDevice.TEMPLATE_PREVIEW)
            captureRequestBuilder?.addTarget(previewSurface)

            //captureRequestBuilder?.set(CaptureRequest.CONTROL_AE_MODE, CaptureRequest.CONTROL_AE_MODE_OFF)
            //captureRequestBuilder?.set(CaptureRequest.CONTROL_MODE, CaptureRequest.CONTROL_MODE_OFF)
            //captureRequestBuilder?.set(CaptureRequest.SENSOR_EXPOSURE_TIME, exposure)

            // camera
            if (cameraController is CameraController1) {
                Log.d(TAG, "Using camera")
                //TODO what if these are null?
                val minExp = cameraController?.cameraFeatures?.min_exposure ?: 0
                val maxExp = cameraController?.cameraFeatures?.max_exposure ?: 0
                val expComp = minExp + progress/100.0*(maxExp - minExp)
                Log.d(TAG, "Min ExpComp $minExp Max ExpComp $maxExp")
                cameraController?.exposureCompensation = expComp.toInt()
            }
            // camera2
            else if (cameraController is CameraController2) {
                Log.d(TAG, "Using camera2")
                /**
                 * CameraController will not allow auto exposure to turn off unless we also set a
                 * manual ISO first
                 */
                cameraController?.setManualISO(true, 100)
                cameraController?.exposureTime = exposure
            }

            Log.d(TAG, "EXP TIME: ${captureRequestBuilder?.get(CaptureRequest.SENSOR_EXPOSURE_TIME)}")
            Log.d(TAG, "Control mode: ${captureRequestBuilder?.get(CaptureRequest.CONTROL_AE_MODE)}")

            captureRequest = captureRequestBuilder?.build()
            Log.d(TAG, "Request: $captureRequest")
            cameraCaptureSession?.setRepeatingRequest(captureRequest, null, backgroundHandler)
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
