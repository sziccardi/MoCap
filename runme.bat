@echo off

REM generate the result pictures

src\Debug\imgpro input\test_sequence_frames\1test_sequence0000001.jpg output\test_sequence_output0000001.jpg -trackMarkers input\reference_image_crop40x40.jpg 1 1

src\Debug\imgpro -svdTests
