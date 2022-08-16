format="TIME: \t%E real,\t%U user,\t%S sys \n CPU: \t%P \t (%Xkb text + %Dkb data: %Mkb max) \n I/O: \t%I inputs + %O outputs \t(%F major + %R minor) pagefaults %W swaps"

log_eval()
{
  cd $1
  echo -e "\nIn $1\n"
  echo "Running: $2"

  if [ "$#" -eq 3 ]; then
    eval "$2" 2>&1 | tee -a $3
    status=${PIPESTATUS[0]}
  elif [ "$#" -eq 2 ]; then
    #/usr/bin/time -f "$format" $2 2>&1
    eval "$2" 2>&1
    status=$?
  else
    echo "Not enough arguments"
    exit 1
  fi

  if [ ! $status == 0 ]; then
    echo "Previous command returned error: $status"
    return 1
  fi
}
