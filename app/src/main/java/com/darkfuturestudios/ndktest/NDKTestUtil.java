package com.darkfuturestudios.ndktest;

import android.app.Activity;
import android.util.Log;
import android.view.Surface;

import com.darkfuturestudios.ndktest.CameraController.CameraController;

public class NDKTestUtil {

    /**
     * If true, CameraController will log
     */
    public static boolean LOG = true;

    private static final String TAG = "NDKTestUtil";

    /**
     * Returns the ROTATION_* enum of the display relative to the natural device orientation.
     */
    private static int getDisplayRotation(Activity activity) {
        // gets the display rotation (as a Surface.ROTATION_* constant), taking into account the getRotatePreviewPreferenceKey() setting
        int rotation = activity.getWindowManager().getDefaultDisplay().getRotation();

        String rotate_preview = "0";
        if (NDKTestUtil.LOG)
            Log.d(TAG, "    rotate_preview = " + rotate_preview);
        if (rotate_preview.equals("180")) {
            switch (rotation) {
                case Surface.ROTATION_0:
                    rotation = Surface.ROTATION_180;
                    break;
                case Surface.ROTATION_90:
                    rotation = Surface.ROTATION_270;
                    break;
                case Surface.ROTATION_180:
                    rotation = Surface.ROTATION_0;
                    break;
                case Surface.ROTATION_270:
                    rotation = Surface.ROTATION_90;
                    break;
                default:
                    break;
            }
        }

        return rotation;
    }

    /**
     * Returns the rotation in degrees of the display relative to the natural device orientation.
     */
    private static int getDisplayRotationDegrees(Activity activity) {
        if (NDKTestUtil.LOG)
            Log.d(TAG, "getDisplayRotationDegrees");
        int rotation = getDisplayRotation(activity);
        int degrees = 0;
        switch (rotation) {
            case Surface.ROTATION_0:
                degrees = 0;
                break;
            case Surface.ROTATION_90:
                degrees = 90;
                break;
            case Surface.ROTATION_180:
                degrees = 180;
                break;
            case Surface.ROTATION_270:
                degrees = 270;
                break;
            default:
                break;
        }
        if (NDKTestUtil.LOG)
            Log.d(TAG, "    degrees = " + degrees);
        return degrees;
    }

    // for the Preview - from http://developer.android.com/reference/android/hardware/Camera.html#setDisplayOrientation(int)
    // note, if orientation is locked to landscape this is only called when setting up the activity, and will always have the same orientation
    public static void setCameraDisplayOrientation(Activity activity, CameraController camera_controller, boolean using_camera2) {
        if (NDKTestUtil.LOG)
            Log.d(TAG, "setCameraDisplayOrientation()");
        if (camera_controller == null) {
            if (NDKTestUtil.LOG)
                Log.d(TAG, "camera not opened!");
            return;
        }
        if (using_camera2) {
            // need to configure the textureview
            // TODO do we need this? Probably not
            // configureTransform();
        } else {
            int degrees = getDisplayRotationDegrees(activity);
            if (NDKTestUtil.LOG)
                Log.d(TAG, "    degrees = " + degrees);
            // note the code to make the rotation relative to the camera sensor is done in camera_controller.setDisplayOrientation()
            camera_controller.setDisplayOrientation(degrees);
        }
    }

}
