# How to create linux wheels

1. Build docker image
2. Login into contrainer
3. Build and repair wheels, manually changing needed python version
 
```bash
sudo docker buildx build .
sudo docker image ls
sudo docker run --name pydegensac2 -v $(pwd):/mnt/pydegensac -it IMG_ID
alias cmake=/usr/local/bin/cmake
/opt/python/cp36-cp36m/bin/python3 setup.py bdist_wheel
auditwheel repair dist/pydegensac-0.1.1-cp36-cp36m-linux_x86_64.whl
```
