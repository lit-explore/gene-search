FROM node:latest

MAINTAINER keith.hughitt@nih.gov

RUN mkdir -p /home/node/app
WORKDIR /home/node/app

# work-around (node v17 x webpack issues)
# https://github.com/webpack/webpack/issues/14532
ENV NODE_OPTIONS="--openssl-legacy-provider"

ADD app/package.json ./
ADD app/package-lock.json ./

RUN npm i

COPY --chown=node:node app/ /home/node/app
