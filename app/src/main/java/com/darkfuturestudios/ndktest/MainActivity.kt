package com.darkfuturestudios.ndktest

import android.content.Context
import android.content.pm.PackageManager
import android.graphics.Bitmap
import android.graphics.SurfaceTexture
import android.hardware.camera2.*
import android.os.Bundle
import android.os.Environment
import android.os.Handler
import android.os.HandlerThread
import android.support.v4.app.ActivityCompat
import android.support.v7.app.AppCompatActivity
import android.util.Log
import android.util.Size
import android.view.Surface
import android.view.TextureView
import android.view.View
import android.widget.Toast
import kotlinx.android.synthetic.main.activity_main.*
import java.io.File
import java.io.FileOutputStream
import java.io.IOException
import java.text.SimpleDateFormat
import java.util.*


class MainActivity : AppCompatActivity() {

    private lateinit var stateCallback: CameraDevice.StateCallback
    private lateinit var textureView: TextureView
    private lateinit var fabTakePhoto: View
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


    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_camera)

        textureView = findViewById(R.id.texture_view)
        fabTakePhoto = findViewById(R.id.fab_take_photo)
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
                openCamera()
            }

            override fun onSurfaceTextureSizeChanged(surfaceTexture: SurfaceTexture, width: Int, height: Int) {

            }

            override fun onSurfaceTextureDestroyed(surfaceTexture: SurfaceTexture): Boolean {
                return false
            }

            override fun onSurfaceTextureUpdated(surfaceTexture: SurfaceTexture) {

            }
        }

        stateCallback = object : CameraDevice.StateCallback() {
            override fun onOpened(cameraDevice: CameraDevice) {
                this@MainActivity.cameraDevice = cameraDevice
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
            setUpCamera()
            openCamera()
        } else {
            textureView.surfaceTextureListener = surfaceTextureListener
        }
    }

    override fun onStop() {
        super.onStop()
        closeCamera()
        closeBackgroundThread()
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
            }
        } catch (e: CameraAccessException) {
            e.printStackTrace()
        }

    }

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

    private fun createPreviewSession() {
        try {
            val surfaceTexture = textureView.surfaceTexture
            surfaceTexture.setDefaultBufferSize(previewSize.width, previewSize.height)
            val previewSurface = Surface(surfaceTexture)
            captureRequestBuilder = cameraDevice?.createCaptureRequest(CameraDevice.TEMPLATE_PREVIEW)
            captureRequestBuilder?.addTarget(previewSurface)

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

        // Used to load the 'native-lib' library on application startup.
        init {
            System.loadLibrary("native-lib")
        }

    }

}