#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"
#include "options.h"

int main(int argc, char ** argv)
{
    Options opt;
    Options_set_defaults(&opt);
    Options_parse_options(argc, argv, &opt);
    Options_print_options(&opt);
    return 0;
}
