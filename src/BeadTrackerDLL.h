//
// Created by Taro Porschke on 6/5/26.
//

#ifndef BEADTRACKERDLL_H
#define BEADTRACKERDLL_H

#ifdef _WIN32
    #ifdef BEADTRACKERDLL_EXPORTS
        #define DLL_API __declspec(dllexport)
    #else
        #define DLL_API __declspec(dllimport)
    #endif
#else
    #define DLL_API
#endif

extern "C" {

    // Called once at experiment start
    // x1,y1 = first bead pixel coords (user clicked in LabVIEW)
    // x2,y2 = second bead pixel coords
    // width, height = frame dimensions
    // returns: opaque handle (pointer cast to int64)
    //          returns 0 if initialization failed
    DLL_API long long createTracker(
        float x1, float y1,
        float x2, float y2,
        int width, int height,
        const char* csv_path,       // Added
        double temperature_K        // Added
    );

    // Called every frame
    // handle = value returned by createTracker
    // pixels = raw grayscale pixel buffer from IMAQ
    // width, height = frame dimensions
    // out_* = pointers LabVIEW provides to receive results
    DLL_API void trackFrame(
        long long handle,
        unsigned char* pixels,
        int width, int height,
        float* out_x1, float* out_y1,
        float* out_x2, float* out_y2,
        float* out_dist_px,
        float* out_dist_um,
        double* out_force_pN,       // Added
        double* out_live_L_um       // Added
    );

    // Called once at experiment end
    DLL_API void destroyTracker(long long handle);

    // Returns processing time of most recent trackFrame call in milliseconds
    DLL_API double getLastFrameTimeMs(long long handle);

    // Resets force calculator buffer for checkpointing
    DLL_API void resetForceCalculator(long long handle);

}
#endif //BEADTRACKERDLL_H