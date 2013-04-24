"""distutils.command.bdist_wininst_amd64

This is a hacked copy of 'bdist_wininst' command to address multiple issues with 
generating python win64 installer on linux 

Should be removed when binutils can generate win64 installer.

"""

__revision__ = "$Id$"

import sys
import os
import string

from distutils.core import Command
from distutils.dir_util import remove_tree
from distutils.errors import DistutilsOptionError, DistutilsPlatformError
from distutils import log
from distutils.util import get_platform

class bdist_wininst_amd64 (Command):

    description = "create an executable installer for MS Windows"

    user_options = [('bdist-dir=', None,
                     "temporary directory for creating the distribution"),
                    ('plat-name=', 'p',
                     "platform name to embed in generated filenames "
                     "(default: %s)" % get_platform()),
                    ('keep-temp', 'k',
                     "keep the pseudo-installation tree around after " +
                     "creating the distribution archive"),
                    ('target-version=', None,
                     "require a specific python version" +
                     " on the target system"),
                    ('no-target-compile', 'c',
                     "do not compile .py to .pyc on the target system"),
                    ('no-target-optimize', 'o',
                     "do not compile .py to .pyo (optimized)"
                     "on the target system"),
                    ('dist-dir=', 'd',
                     "directory to put final built distributions in"),
                    ('bitmap=', 'b',
                     "bitmap to use for the installer instead of python-powered logo"),
                    ('title=', 't',
                     "title to display on the installer background instead of default"),
                    ('skip-build', None,
                     "skip rebuilding everything (for testing/debugging)"),
                    ('install-script=', None,
                     "basename of installation script to be run after"
                     "installation or before deinstallation"),
                    ('pre-install-script=', None,
                     "Fully qualified filename of a script to be run before "
                     "any files are installed.  This script need not be in the "
                     "distribution"),
                    ('user-access-control=', None,
                     "specify Vista's UAC handling - 'none'/default=no "
                     "handling, 'auto'=use UAC if target Python installed for "
                     "all users, 'force'=always use UAC"),
                   ]

    boolean_options = ['keep-temp', 'no-target-compile', 'no-target-optimize',
                       'skip-build']

    def initialize_options (self):
        self.bdist_dir = None
        self.plat_name = None
        self.keep_temp = 0
        self.no_target_compile = 0
        self.no_target_optimize = 0
        self.target_version = None
        self.dist_dir = None
        self.bitmap = None
        self.title = None
        self.skip_build = None
        self.install_script = None
        self.pre_install_script = None
        self.user_access_control = None

    # initialize_options()


    def finalize_options (self):
        self.set_undefined_options('bdist', ('skip_build', 'skip_build'))
        self.plat_name = 'win-x86_64'

        if self.bdist_dir is None:
            if self.skip_build and self.plat_name:
                # If build is skipped and plat_name is overridden, bdist will
                # not see the correct 'plat_name' - so set that up manually.
                bdist = self.distribution.get_command_obj('bdist')
                bdist.plat_name = self.plat_name
                # next the command will be initialized using that name
            bdist_base = self.get_finalized_command('bdist').bdist_base
            self.bdist_dir = os.path.join(bdist_base, 'wininst')

        if not self.target_version:
            self.target_version = ""

        self.set_undefined_options('bdist',
                                   ('dist_dir', 'dist_dir'),
                                   ('plat_name', 'plat_name'),
                                  )

        if self.install_script:
            for script in self.distribution.scripts:
                if self.install_script == os.path.basename(script):
                    break
            else:
                raise DistutilsOptionError, \
                      "install_script '%s' not found in scripts" % \
                      self.install_script
    # finalize_options()


    def run (self):
        if (sys.platform != "win32" and
            (self.distribution.has_ext_modules() or
             self.distribution.has_c_libraries())):
            raise DistutilsPlatformError \
                  ("distribution contains extensions and/or C libraries; "
                   "must be compiled on a Windows 32 platform")

        if not self.skip_build:
            self.run_command('build')

        install = self.reinitialize_command('install', reinit_subcommands=1)
        install.root = self.bdist_dir
        install.skip_build = self.skip_build
        install.warn_dir = 0
        install.plat_name = self.plat_name

        install_lib = self.reinitialize_command('install_lib')
        # we do not want to include pyc or pyo files
        install_lib.compile = 0
        install_lib.optimize = 0

        if self.distribution.has_ext_modules():
            # If we are building an installer for a Python version other
            # than the one we are currently running, then we need to ensure
            # our build_lib reflects the other Python version rather than ours.
            # Note that for target_version!=sys.version, we must have skipped the
            # build step, so there is no issue with enforcing the build of this
            # version.
            target_version = self.target_version
            if not target_version:
                assert self.skip_build, "Should have already checked this"
                target_version = sys.version[0:3]
            plat_specifier = ".%s-%s" % (self.plat_name, target_version)
            build = self.get_finalized_command('build')
            build.build_lib = os.path.join(build.build_base,
                                           'lib' + plat_specifier)

        # Use a custom scheme for the zip-file, because we have to decide
        # at installation time which scheme to use.
        for key in ('purelib', 'platlib', 'headers', 'scripts', 'data'):
            value = str.upper(key)
            if key == 'headers':
                value = value + '/Include/$dist_name'
            setattr(install,
                    'install_' + key,
                    value)

        log.info("installing to %s", self.bdist_dir)
        install.ensure_finalized()

        # avoid warning of 'install_lib' about installing
        # into a directory not in sys.path
        sys.path.insert(0, os.path.join(self.bdist_dir, 'PURELIB'))

        install.run()

        del sys.path[0]

        # And make an archive relative to the root of the
        # pseudo-installation tree.
        from tempfile import mktemp
        archive_basename = mktemp()
        fullname = self.distribution.get_fullname()
        arcname = self.make_archive(archive_basename, "zip",
                                    root_dir=self.bdist_dir)
        # create an exe containing the zip-file
        self.create_exe(arcname, fullname, self.bitmap)
        pyversion = 'any'
        self.distribution.dist_files.append(('bdist_wininst', pyversion,
                                             self.get_installer_filename(fullname)))
        # remove the zip-file again
        log.debug("removing temporary file '%s'", arcname)
        os.remove(arcname)

        if not self.keep_temp:
            remove_tree(self.bdist_dir, dry_run=self.dry_run)

    # run()

    def get_inidata (self):
        # Return data describing the installation.

        lines = []
        metadata = self.distribution.metadata

        # Write the [metadata] section.
        lines.append("[metadata]")

        # 'info' will be displayed in the installer's dialog box,
        # describing the items to be installed.
        info = (metadata.long_description or '') + '\n'

        # Escape newline characters
        def escape(s):
            return str.replace(s, "\n", "\\n")

        for name in ["author", "author_email", "description", "maintainer",
                     "maintainer_email", "name", "url", "version"]:
            data = getattr(metadata, name, "")
            if data:
                info = info + ("\n    %s: %s" % \
                               (str.capitalize(name), escape(data)))
                lines.append("%s=%s" % (name, escape(data)))

        # The [setup] section contains entries controlling
        # the installer runtime.
        lines.append("\n[Setup]")
        if self.install_script:
            lines.append("install_script=%s" % self.install_script)
        lines.append("info=%s" % escape(info))
        lines.append("target_compile=%d" % (not self.no_target_compile))
        lines.append("target_optimize=%d" % (not self.no_target_optimize))
        if self.target_version:
            lines.append("target_version=%s" % self.target_version)
        if self.user_access_control:
            lines.append("user_access_control=%s" % self.user_access_control)

        title = self.title or self.distribution.get_fullname()
        lines.append("title=%s" % escape(title))
        import time
        import distutils
        build_info = "Built %s with distutils-%s" % \
                     (time.ctime(time.time()), distutils.__version__)
        lines.append("build_info=%s" % build_info)
        return "\n".join(lines)

    # get_inidata()

    def create_exe (self, arcname, fullname, bitmap=None):
        import struct

        self.mkpath(self.dist_dir)

        cfgdata = self.get_inidata()

        installer_name = self.get_installer_filename(fullname)
        self.announce("creating %s" % installer_name)

        if bitmap:
            bitmapdata = open(bitmap, "rb").read()
            bitmaplen = len(bitmapdata)
        else:
            bitmaplen = 0

        file = open(installer_name, "wb")
        file.write(self.get_exe_bytes())
        if bitmap:
            file.write(bitmapdata)

        # Convert cfgdata from unicode to ascii, mbcs encoded
        try:
            unicode
        except NameError:
            pass
        else:
            if isinstance(cfgdata, unicode):
                cfgdata = cfgdata.encode("ascii")

        # Append the pre-install script
        cfgdata = cfgdata + "\0".encode("ascii")
        if self.pre_install_script:
            script_data = open(self.pre_install_script, "r").read()
            cfgdata = cfgdata + script_data + "\n\0".encode("ascii")
        else:
            # empty pre-install script
            cfgdata = cfgdata + "\0".encode("ascii")
        file.write(cfgdata)

        # The 'magic number' 0x1234567B is used to make sure that the
        # binary layout of 'cfgdata' is what the wininst.exe binary
        # expects.  If the layout changes, increment that number, make
        # the corresponding changes to the wininst.exe sources, and
        # recompile them.
        header = struct.pack("<iii",
                             0x1234567B,       # tag
                             len(cfgdata),     # length
                             bitmaplen,        # number of bytes in bitmap
                             )
        file.write(header)
        file.write(open(arcname, "rb").read())

    # create_exe()

    def get_installer_filename(self, fullname):
        # Factored out to allow overriding in subclasses
        if self.target_version:
            # if we create an installer for a specific python version,
            # it's better to include this in the name
            installer_name = os.path.join(self.dist_dir,
                                          "%s.%s-py%s.exe" %
                                           (fullname, self.plat_name, self.target_version))
        else:
            installer_name = os.path.join(self.dist_dir,
                                          "%s.%s.exe" % (fullname, self.plat_name))
        return installer_name
    # get_installer_filename()

    def get_exe_bytes (self):
        directory = os.path.dirname(__file__)
        filename = os.path.join(directory, "tools", "wininst-9.0-amd64.exe")
        f = open(filename, "rb")
        try:
            return f.read()
        finally:
            f.close()
# class bdist_wininst
