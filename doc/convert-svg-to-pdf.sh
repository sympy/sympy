#!/bin/bash
#
# Convert an SVG file to a PDF file by using headless Chrome.
#
# This is based on
# https://gist.github.com/s417-lama/84bf66de1096c4587e8187092fb41684. See also
# https://superuser.com/questions/381125/how-do-i-convert-an-svg-to-a-pdf-on-linux
# and
# https://stackoverflow.com/questions/44970113/how-can-i-change-paper-size-in-headless-chrome-print-to-pdf
#
# Note: if you end up with a pdf with a blank page, you need to edit the SVG
# file to trim down the width and/or height.

if [ $# -ne 2 ]; then
  echo "Usage: ./convert-svg-to-pdf.sh input.svg output.pdf" 1>&2
  exit 1
fi

INPUT=$1
OUTPUT=$2

HTML="
<html>
  <head>
    <style>
body {
  margin: 0;
}
    </style>
    <script>
function init() {
  const element = document.getElementById('targetsvg');
  const positionInfo = element.getBoundingClientRect();
  const height = positionInfo.height;
  const width = positionInfo.width;
  const style = document.createElement('style');
  style.innerHTML = \`@page {margin: 0; size: \${width}px \${height}px}\`;
  document.head.appendChild(style);
}
window.onload = init;
    </script>
  </head>
  <body>
    <img id=\"targetsvg\" src=\"${INPUT}\">
  </body>
</html>
"

tmpfile=$(mktemp XXXXXX.html)
trap "rm -f $tmpfile" EXIT
echo $HTML > "$tmpfile"

if command -v google-chrome &> /dev/null
then
    CHROME=google-chrome
elif command -v "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" &> /dev/null
then
    CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
elif command -v chrome &> /dev/null
then
    CHROME=chrome
else
    echo "Could not find a valid command for Google Chrome. Make sure Chrome is installed."
    exit 1
fi

"$CHROME" --headless --disable-gpu --print-to-pdf="$OUTPUT" "$tmpfile"
