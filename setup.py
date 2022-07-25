
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.readlines()


setup(
        name ='B2RIO',
        version ='0.0.1',
        author ='Gaston E. Zanitti',
        author_email ='gzanitti@gmail.com',
        url ='https://github.com/gzanitti/B2RIO',
        description ='Brain t(w)o Reverse Inference Ontology',
        license ='MIT',
        packages = find_packages(),
        python_requires='>3.8',
        entry_points ={
            'console_scripts': [
                'b2rio = b2rio.b2rio:run',
                'b2rio_prob = b2rio.b2rio:run_probabilistic'
            ]
        },
        classifiers =[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        install_requires = requirements,
        zip_safe = False
)