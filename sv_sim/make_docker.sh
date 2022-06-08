################################################################################
## build the micronap docker
################################################################################

# This is a very simpel thing for building a docker
# add the tag version as input

echo $PWD

#log_eval $PWD "docker image rm -f micronap:$1"
VERSION=$1

docker build --build-arg VERSION=$VERSION . --tag sv_sim:$VERSION
docker tag sv_sim:$VERSION vrohnie/sv_sim:$VERSION
# docker push vrohnie/sv_sim:$VERSION
