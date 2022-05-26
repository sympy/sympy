from __future__ import annotations
"""
    Converts SVG images to PDF using chrome in case the builder does not
    support SVG images natively (e.g. LaTeX).

"""
import subprocess

from sphinx.errors import ExtensionError
from sphinx.locale import __
from sphinx.transforms.post_transforms.images import ImageConverter
from sphinx.util import logging
from errno import ENOENT, EPIPE, EINVAL
import os
import platform
import tempfile
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

    def chrome_command(self) -> str:
        if platform.win32_ver()[0]:
            return os.path.join(os.environ['PROGRAMW6432'],
             "Google\Chrome\Application\chrome.exe")
        elif platform.mac_ver()[0]:
            return "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        elif platform.libc_ver()[0]:
            return "google-chrome"
        else:
            logger.error("Not able to find suitable chromium or chrome command for this os")

    def chromium_command(self) -> str:
        if platform.win32_ver()[0]:
            return os.path.join(os.environ['PROGRAMW6432'],
             "chrome-win/chrome.exe")
        elif platform.mac_ver()[0]:
            return "/Applications/Chromium.app/Contents/MacOS/Chromium"
        elif platform.libc_ver()[0]:
            return "chromium-browser"

    def command_runner(self, chrome, _to, temp_name):
        args = [chrome, '--headless', '--disable-gpu', '--disable-software-rasterizer', '--print-to-pdf=' + _to, temp_name]
        logger.debug('Invoking %r ...', args)
        p = subprocess.run(
            args, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        return p

    def convert(self, _from: str, _to: str) -> bool:
        """Converts the image from SVG to PDF using chrome."""

        HTML = "<html ><head ><style >body {margin: 0; }</style ><script >function init() {const element = document.getElementById('targetsvg');const positionInfo = element.getBoundingClientRect();const height = positionInfo.height;const width = positionInfo.width;const style = document.createElement('style');style.innerHTML = `@page {margin: 0; size: ${width}px ${height}px}`;document.head.appendChild(style); }window.onload = init;</script ></head><body><img id=\"targetsvg\" src=\"%s\"></body></html>" % (_from)
        try:
            temp = tempfile.NamedTemporaryFile(suffix=".html", mode='w', delete=False)
            with open(temp.name, 'w+') as f:
                f.write(HTML)
                f.seek(0)
            chromium = self.chromium_command()
            p = self.command_runner(chromium, _to, temp.name)

        except OSError as err:
            if err.errno != ENOENT:  # No such file or directory
                raise
            try:
                chrome = self.chrome_command()
                p = self.command_runner(chrome, _to, temp.name)
            except OSError as err:
                logger.error('converter command cannot be run. chrome is not a internal or external command')
                return False

        try:
            stdout, stderr = p.stdout, p.stderr
            temp.close()
            os.unlink(temp.name)
        except (OSError, IOError) as err:
            if err.errno not in (EPIPE, EINVAL):
                raise
            stdout, stderr = p.stdout, p.stderr
        if p.returncode != 0:
            raise ExtensionError(__('converter exited with error:\n'
                                    '[stderr]\n%s\n[stdout]\n%s') %
                                 (stderr, stdout))

        return True


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_post_transform(Converter)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
