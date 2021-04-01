from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='rxn_network',
        version='0.0.1',
        description='Stripped down reaction network code',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/GENESIS-EFRC/reaction-network',
        author=['Matthew McDermott'],
        author_email=['mcdermott@lbl.gov'],
        license='MIT',
        packages=['rxn_network'],
        zip_safe=False,
        install_requires=[],
        extras_require={},
        classifiers=[],
        test_suite='',
        tests_require=[],
        scripts=[]
    )
