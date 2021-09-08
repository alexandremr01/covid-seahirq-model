FROM python:3.7-buster
RUN apt-get -y update && apt-get -y install lua5.1 lua-socket lua-sec premake4

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt
ENV PYTHONPATH /app/bin

ENTRYPOINT bash