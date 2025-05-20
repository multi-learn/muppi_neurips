# -*- coding: utf-8 -*-

# Extracting requirements from requirements.txt
import os, re
import io
# from Cython.Build import cythonize
from setuptools import setup, find_packages
from distutils.command.clean import clean as _clean
from distutils.dir_util import remove_tree
from distutils.command.sdist import sdist
from setuptools import setup, find_packages


USE_COPYRIGHT = False
try:
    from copyright import writeStamp, eraseStamp
except ImportError:
    USE_COPYRIGHT = False

# --------------------------------------------------------------------
# Clean target redefinition - force clean everything supprimer de la liste '^core\.*$',
relist = ['^.*~$', '^#.*#$', '^.*\.aux$', '^.*\.pyc$', '^.*\.o$']
reclean = []
###################
# Get MUPPi version
####################
def get_version():
    v_text = open('VERSION').read().strip()
    v_text_formted = '{"' + v_text.replace('\n', '","').replace(':', '":"')
    v_text_formted += '"}'
    v_dict = eval(v_text_formted)
    return v_dict["muppi"]

########################
# Set MuPPI __version__
########################
def set_version(muppi_dir, version):
    filename = os.path.join(muppi_dir, '__init__.py')
    buf = ""
    for line in open(filename, "rb"):
        if not line.decode("utf8").startswith("__version__ ="):
            buf += line.decode("utf8")
    f = open(filename, "wb")
    f.write(buf.encode("utf8"))
    f.write(('__version__ = "%s"\n' % version).encode("utf8"))

for restring in relist:
    reclean.append(re.compile(restring))


def wselect(args, dirname, names):
    for n in names:
        for rev in reclean:
            if (rev.match(n)):
                os.remove("%s/%s" %(dirname, n))
        break
class clean(_clean):
    def walkAndClean(self):
        os.walk("..", wselect, [])
        pass

    def run(self):
        clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('multimodal'):
            for filename in filenames:
                if (filename.endswith('.so') or
                        filename.endswith('.pyd') or
                        filename.endswith('.dll') or
                        filename.endswith('.pyc')):
                    os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == '__pycache__':
                    shutil.rmtree(os.path.join(dirpath, dirname))

class m_sdist(sdist):
    """ Build source package

    WARNING : The stamping must be done on an default utf8 machine !
    """
    def run(self):
        if USE_COPYRIGHT:
            writeStamp()
            sdist.run(self)
            # eraseStamp()
        else:
            sdist.run(self)

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


##########################
# File path read command
##########################
def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with io.open(os.path.join(*paths), 'r', encoding='utf-8') as f:
        return f.read()


def setup_package():
    group = 'multi-learn'
    name='MuPPI_dataset'
    version = get_version()
    muppi_dir = ''
    long_description=read('README.md')
    packages = find_packages()
    setup(version=version,
        packages=packages,
        long_description=long_description,
        include_package_data=True,
        license="BSD-3-Clause",
        license_files="LICENSE"
    )


if __name__ == "__main__":
    setup_package()
