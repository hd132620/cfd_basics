#include "setting.h"
#include "navier_stokes.h"

int main(void) {
    ns_init();
    ns_solve();
    ns_finalize();
    return 0;
}