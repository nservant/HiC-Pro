
if [ $# != 1 ]; then echo "usage: $0 CONFIG"; exit 1; fi

dir=$(dirname $0)

. $dir/hic.inc.sh

config=$1
read_config $config

env
