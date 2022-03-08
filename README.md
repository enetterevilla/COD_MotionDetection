# COD_MotionDetection
COD_Motion_Detection_ForGitHub.m is the adaptive data-driven motion detection algorithm developed by our group at the Yale PET Center. This is the motion detection used in the paper "Adaptive Data-driven Motion Detection and Optimized Correction for Brain PET" submitted to NeuroImage. This article has been accepted for publication and undergone full peer review. Please cite this article as https://doi.org/10.1016/j.neuroimage.2022.119031. 

# Tutorial
The COD file "COD_sample.cod" contains the averaged location of the events every 1 second. This algorithm automatically detects motion based on the COD trace and is adaptive to different tracers and noise levels without user-defined parameters. 

1. Download the .m and .cod files.
2. You can change the start and end time, alongside with alpha, n_max and scout_segments. Note: The recommended n_max for a 90-min FDG scan is 300. 

# Results
The most important result would be "MFF_union_all.mat". It contains the motion-free frame (MFF) information (start and end time) to be used for frame reconstruction and registration later on.
