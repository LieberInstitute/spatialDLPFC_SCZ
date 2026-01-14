To install a local python environment with tensorQTL and rpy2, run the following commands in your Linux terminal (assuming CUDA drivers are installed when using GPU acceleration, and uv present):

```bash
uv venv --python 3.12
source ./venv/bin/activate
uv pip install -U pip setuptools wheel
## check what CUDA version is supported on your system: https://pytorch.org/get-started/locally/
## example for CUDA 12.6:
uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cu126
uv pip install tensorqtl Pgenlib rpy2 pandas scipy numpy pandas_plink
```

For vscode/positron users, an interactive kernel/console can be also installed in the same environment:

```bash
uv pip install ipykernel ipython
```
