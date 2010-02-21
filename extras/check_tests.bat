@echo off

rem This file is used to test IT++ library instalation in Windows
rem It is assumed that the current folder is extras, the *.exe files are in the "bin" folder
rem and the reference files are present in "..\tests\*.ref" folder.
rem The results of the tests are written into "check_tests.txt" file.

rem generate results from test executables
cd ..\win32\bin
for /f %%a in ('dir /b *.exe') do echo Running %%~na.exe & %%~na.exe > %%~na.tmp

rem compare results with reference files
cd ..\..\extras
fc /L ..\win32\bin\*.tmp ..\tests\*.ref > check_test.txt

rem delete results files
del ..\win32\bin\*.tmp

echo Open "check_test.txt" file to verify tests results

pause
