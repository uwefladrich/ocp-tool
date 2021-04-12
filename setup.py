import setuptools


setuptools.setup(
    name='ocp-tool',
    author='Jan Streffing',
    author_email='jan.streffing@awi.de',
    description='Tool to generate OASIS files for coupling OpenIFS, FESOM2, and NEMO',
    url='https://github.com/JanStreffing/ocp-tool',
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    entry_points={
        'scriptengine.tasks': [
            'ocpt.main = ocp_tool.scriptengine_task:OCPTool',
        ],
    },
)
