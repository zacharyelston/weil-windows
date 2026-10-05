# Pinned environment

```bash
docker build -t weil-windows-env:py3.14.4 -f docker/Dockerfile docker
docker/run.sh scripts/decay_claims_check.py
```

`run.sh` mounts the checkout at `/rh2` and runs Python inside the image. The image installs `requirements.txt`: mpmath 1.4.1, gmpy2 2.3.1, python-flint 0.9.0 and cypari2 2.2.4. The last is needed only by `scripts/borromean`.
