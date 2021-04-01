from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='y2mn2o7_selectivity',
        version='0.0.1',
        description='Code/exp_data for thermochemical connectivity manuscript.',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/GENESIS-EFRC/thermochemical-connectivity',
        author=['Matthew McDermott', 'Paul Todd'],
        author_email=['mcdermott@lbl.gov', 'ptodd@nrel.gov'],
        license='MIT',
        packages=['thermochemical_connectivity'],
        zip_safe=False,
        install_requires=[],
        extras_require={},
        classifiers=[],
        test_suite='',
        tests_require=[],
        scripts=[]
    )
