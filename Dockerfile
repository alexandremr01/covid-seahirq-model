FROM python:3.7-buster
RUN apt-get -y update && apt-get -y install cmake

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt
ENV PYTHONPATH /app/build
ENV EIGEN3_INCLUDE_DIR /app/3rdparty/eigen
ENV CMAKE_PREFIX_PATH /usr/local/lib/python3.7/site-packages/pybind11/share/cmake/pybind11

ENTRYPOINT bash