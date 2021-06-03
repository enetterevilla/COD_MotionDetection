# COD_MotionDetection
COD_Motion_Detection_ForGitHub.m is the adaptive data-driven motion detection algorithm developed by our group at the Yale PET Center. This is the motion detection used in the paper "Adaptive Data-driven Motion Detection and Optimized Correction for Brain PET" submitted to NeuroImage. The COD file "COD_sample.cod" contains the averaged location of the events every 1 second. This algorithm automatically detects motion based on the COD trace and is adaptive to different tracers and noise levels without user-defined parameters. 

# Tutorial
1. Download the .m and .cod files.
2. You can change the start and end time, alongside with alpha, n_max and scout_segments. Note: The recommended n_max for a 90-min FDG scan is 300. 
