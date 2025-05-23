#!/bin/sh
set -ex
if [ ${target_platform} == "linux-ppc64le" ]; then
  exit 0
fi
if [ "${JAVA_HOME}" != "${PREFIX}/lib/jvm" ] && [ "${JAVA_HOME}" != "${PREFIX}/Library/lib/jvm" ]; then
    echo "ERROR: JAVA_HOME (${JAVA_HOME}) not equal to PREFIX (${PREFIX}/lib/jvm or ${PREFIX}/Library/lib/jvm)"
  exit 1
fi

pushd test-nio
  javac TestFilePaths.java
  jar cfm TestFilePaths.jar manifest.mf TestFilePaths.class
  java -jar TestFilePaths.jar TestFilePaths.java
popd


set -ex



java -version
${JAVA_HOME}/bin/java -version
./test-jni.sh
test -f ${JAVA_HOME}/release
exit 0
