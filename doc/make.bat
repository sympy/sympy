@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SPHINXOPTS=-W --keep-going
set SOURCEDIR=src
set BUILDDIR=_build
set PORT=3000
set HOST=localhost
set ALLSPHINXOPTSlatex=-d %BUILDDIR%\doctrees-latex -D latex_element.papersize=%PAPER%paper %SPHINXOPTS% %SOURCEDIR%
SETLOCAL
set LATEXMKOPTS=-halt-on-error -xelatex

if "%1" == "html" goto html
if "%1" == "help" goto help
if "%1" == "livehtml" goto livehtml
if "%1" == "pdf" goto pdf

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.https://www.sphinx-doc.org/
	exit /b 1
)

:html
if NOT exist %BUILDDIR%/ (
	mkdir %BUILDDIR%\logo
	python ./generate_logos.py -d
)
%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:livehtml
if NOT exist %BUILDDIR%/ (
	mkdir %BUILDDIR%\logo
	python ./generate_logos.py -d
)
sphinx-autobuild --port %PORT% --host %HOST% --open-browser --watch .. %SOURCEDIR% %BUILDDIR%\html
goto end

:pdf
if NOT exist %BUILDDIR%/logo (
	mkdir %BUILDDIR%\logo %BUILDDIR%\latex %BUILDDIR%\doctrees
	python ./generate_logos.py -d
)
%SPHINXBUILD% -b latex %ALLSPHINXOPTSlatex% %BUILDDIR%\latex
pushd %BUILDDIR%\latex
latexmk -xelatex -silent -f
popd
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
echo   [94mlivehost[0m    to create and host HTML on [94mhttp://localhost:3000[0m

:end
popd
