from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='y2mn2o7_selectivity',
        version='0.0.1',
        description='Code and exp_data for Y2Mn2O7 selectivity manuscript.',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        url='https://github.com/GENESIS-EFRC/y2mn2o7-selectivity',
        author=['Matthew McDermott', 'Paul Todd'],
        author_email=['mcdermott@lbl.gov', 'ptodd@nrel.gov'],
        license='MIT',
        packages=['y2mn2o7_selectivity'],
        zip_safe=False,
        install_requires=[],
        extras_require={},
        classifiers=[],
        test_suite='',
        tests_require=[],
        scripts=[]
    )
