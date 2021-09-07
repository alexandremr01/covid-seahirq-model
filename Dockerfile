FROM python:3.7-buster
RUN apt-get -y update && apt-get -y install lua5.1 lua-socket lua-sec

WORKDIR /usr/local/bin
RUN curl -s -L https://github.com/premake/premake-core/releases/download/v5.0.0-alpha16/premake-5.0.0-alpha16-linux.tar.gz | tar -xzf - && \
    chmod +x premake5

WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

ENTRYPOINT bash