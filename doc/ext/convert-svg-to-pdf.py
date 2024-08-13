from __future__ import annotations
"""
    Converts SVG images to PDF using chrome in case the builder does not
    support SVG images natively (e.g. LaTeX).

"""

from sphinx.transforms.post_transforms.images import ImageConverter
from sphinx.util import logging
import os
import platform
from typing import Any  # NOQA
from sphinx.application import Sphinx  # NOQA


logger = logging.getLogger(__name__)


class Converter(ImageConverter):
    conversion_rules = [
        ('image/svg+xml', 'application/pdf'),
    ]

    def is_available(self) -> bool:
        """Confirms if converter is available or not."""
        return True

    def chrome_command(self) -> str | None:
        if platform.win32_ver()[0]:
            if os.system("where chrome") == 0:
                return "chrome"
            path = os.path.join(os.environ["PROGRAMW6432"], "Google\\Chrome\\Application\\chrome.exe")
            if os.path.exists(path):
                return f'"{path}"'
            return None
        if os.system("chrome --version") == 0:
            return "chrome"
        if platform.mac_ver()[0]:
            return "'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'"
        elif platform.libc_ver()[0]:
            return "google-chrome"
        return None

    def chromium_command(self) -> str | None:
        if platform.win32_ver()[0]:
            if os.system("where chromium") == 0:
                return "chromium"
            path = os.path.join(os.environ["PROGRAMW6432"], "Chromium\\Application\\chrome.exe")
            if os.path.exists(path):
                return f'"{path}"'
            return None
        if os.system("chromium --version") == 0:
            return "chromium"
        if platform.mac_ver()[0]:
            path = "/Applications/Chromium.app/Contents/MacOS/Chromium"
            if os.path.exists(path):
                return path
        elif platform.libc_ver()[0]:
            if os.system("chromium-browser --version") == 0:
                return "chromium-browser"
        return None


    def command_runner(self, chrome: str | None, _to: str, temp_name: str) -> int:
        if not chrome:
            return 1
        command = f'{chrome} --headless --disable-gpu --disable-software-rasterizer --print-to-pdf={_to} {temp_name}'
        logger.info(command)
        return os.system(command)

    def convert(self, _from: str, _to: str) -> bool:
        """Converts the image from SVG to PDF using chrome."""
        with open(_from, 'r') as f:
            svg = f.read()

        HTML = "<html><head><style>body {margin: 0; }</style><script>function init() {const element = document.querySelector('svg');const positionInfo = element.getBoundingClientRect();const height = positionInfo.height;const width = positionInfo.width;const style = document.createElement('style');style.innerHTML = `@page {margin: 0; size: ${width}px ${height}px}`;document.head.appendChild(style); }window.onload = init;</script></head><body>%s</body></html>" % (svg)
        temp_name = f'{_from}.html'
        with open(temp_name, 'w') as f:
            f.write(HTML)

        chromium = self.chromium_command()
        code = self.command_runner(chromium, _to, temp_name)
        if code != 0:
            chrome = self.chrome_command()
            code = self.command_runner(chrome, _to, temp_name)
        if code != 0:
            logger.error('Fail to convert svg to pdf. Make sure Chromium or Chrome is installed.')
            exit(1)
        return True


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_post_transform(Converter)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
