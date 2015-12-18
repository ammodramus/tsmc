#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"
#include "options.h"

int main(int argc, char ** argv)
{
    Options opt;
    Options_set_defaults(&opt);
    Options_parse_options(&opt, argc, argv);
    Options_print_options(&opt);
    return 0;
}
