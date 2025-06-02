function log ()
{
   printf "\n\n"
   printf "\e[93m--~*~--\e[m "
   printf "\e[30;103m $@ \e[m"
   printf " \e[93m--~*~--\e[m "
   printf "\n"
}

function build_mmpio ()
{
    go build -C ../../src -v -o ../mmpio
}
