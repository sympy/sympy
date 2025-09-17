@ECHO OFF
SETLOCAL ENABLEEXTENSIONS ENABLEDELAYEDEXPANSION

REM ==== Base config ==========================================================
IF "%SPHINXBUILD%"=="" SET "SPHINXBUILD=sphinx-build"
REM Default: do NOT fail build on warnings (use `html-strict` to enforce -W)
SET "SPHINXOPTS="
SET "STRICTOPTS=-W --keep-going"
SET "SOURCEDIR=src"
SET "BUILDDIR=_build"
SET "PORT=3000"
SET "HOST=localhost"
SET "LATEXMKOPTS=-halt-on-error -xelatex"

REM Precompute latex opts (no -W here either)
SET "ALLSPHINXOPTS_LATEX=-d \"%BUILDDIR%\doctrees-latex\" %SPHINXOPTS% \"%SOURCEDIR%\""

REM Pick Python launcher
WHERE py >NUL 2>&1
IF ERRORLEVEL 1 (
  SET "PY=python"
) ELSE (
  SET "PY=py -3"
)

REM ==== Targets ==============================================================
IF "%~1"=="" GOTO help
IF /I "%~1"=="help" GOTO help
IF /I "%~1"=="html" GOTO html
IF /I "%~1"=="html-strict" GOTO html_strict
IF /I "%~1"=="livehtml" GOTO livehtml
IF /I "%~1"=="pdf" GOTO pdf

ECHO Unknown target "%~1". Run "make.bat help".
EXIT /B 2

REM ---- utility: check sphinx is available ----
:check_sphinx
%SPHINXBUILD% >NUL 2>NUL
IF ERRORLEVEL 9009 (
  ECHO.
  ECHO The 'sphinx-build' command was not found. Make sure you have Sphinx
  ECHO installed, then set the SPHINXBUILD environment variable to point
  ECHO to the full path of the 'sphinx-build' executable. Alternatively you
  ECHO may add the Sphinx directory to PATH.
  ECHO.
  ECHO If you don't have Sphinx installed, grab it from:
  ECHO https://www.sphinx-doc.org/
  ECHO.
  EXIT /B 1
)
EXIT /B 0

REM ---- utility: ensure build dir & logos exist ----
:prepare_builddir
IF NOT EXIST "%BUILDDIR%\logo" (
  MKDIR "%BUILDDIR%\logo"
  %PY% "./generate_logos.py" -d
)
EXIT /B 0

REM ---- html (warnings do NOT fail) ----
:html
CALL :check_sphinx || EXIT /B 1
CALL :prepare_builddir
%SPHINXBUILD% -M html "%SOURCEDIR%" "%BUILDDIR%" %SPHINXOPTS% %O%
GOTO end

REM ---- html-strict (warnings FAIL the build) ----
:html_strict
CALL :check_sphinx || EXIT /B 1
CALL :prepare_builddir
%SPHINXBUILD% -M html "%SOURCEDIR%" "%BUILDDIR%" %STRICTOPTS% %SPHINXOPTS% %O%
GOTO end

REM ---- livehtml dev server ----
:livehtml
CALL :check_sphinx || EXIT /B 1
CALL :prepare_builddir
sphinx-autobuild --port "%PORT%" --host "%HOST%" --open-browser --watch ".." "%SOURCEDIR%" "%BUILDDIR%\html"
GOTO end

REM ---- pdf via LaTeX (warnings do NOT fail) ----
:pdf
CALL :check_sphinx || EXIT /B 1
IF NOT EXIST "%BUILDDIR%\logo" (
  MKDIR "%BUILDDIR%\logo" "%BUILDDIR%\latex" "%BUILDDIR%\doctrees"
  %PY% "./generate_logos.py" -d
)
%SPHINXBUILD% -b latex %ALLSPHINXOPTS_LATEX% "%BUILDDIR%\latex"
PUSHD "%BUILDDIR%\latex"
latexmk %LATEXMKOPTS% -silent -f
POPD
GOTO end

REM ---- help ----
:help
ECHO Usage: make.bat [TARGET]
ECHO.
ECHO Targets:
ECHO   html           Build Sphinx HTML (warnings do NOT fail the build)
ECHO   html-strict    Build HTML and *fail on warnings* (-W --keep-going)
ECHO   livehtml       Auto-reload docs server on http://%HOST%:%PORT%
ECHO   pdf            Build PDF via LaTeX
ECHO.
EXIT /B 0

:end
ENDLOCAL
EXIT /B %ERRORLEVEL%
