# native numerical environment: environment only, code bind-mounted at /work.
# Build:  docker build -t weil-windows-env:py3.14.4 -f docker/Dockerfile docker
# Run:    docker/run.sh scripts/grid_certificate_flint.py --mode zhu --a 0.8 --tsharp 200 --nmodes 120
FROM python:3.14.4-slim-bookworm
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt \
 && python -c "import mpmath, gmpy2, flint; assert mpmath.libmp.BACKEND == 'gmpy'; print('mpmath', mpmath.__version__, 'gmpy2', gmpy2.version(), 'flint', flint.__version__)"
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /work
ENTRYPOINT ["python"]
