For /R HRLQ_outputs\master %%? IN (*.txt) do (
echo %%~n?%%~x? >> %%~dp?outputs\test_outputs.mxt
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o %%~dp?outputs\%%~n?_output.mxt >> %%~dp?outputs\test_outputs.mxt
)
forfiles /S /M *.mxt /C "cmd /c rename @file @fname.txt"