For /R test_instances %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o %%~dp?%%~n?_output.mxt
)
forfiles /S /M *.mxt /C "cmd /c rename @file @fname.txt"