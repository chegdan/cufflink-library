#include <thrust/version.h>
#include <cusp/version.h>
#include "version.h"
#include <iostream>
        
int main(void)
{
    int thrust_major = THRUST_MAJOR_VERSION;
    int thrust_minor = THRUST_MINOR_VERSION;
        
    int cusp_major = CUSP_MAJOR_VERSION;
    int cusp_minor = CUSP_MINOR_VERSION;
    int cusp_subminor = CUSP_SUBMINOR_VERSION;

    int cufflink_major = CUFFLINK_MAJOR_VERSION;
    int cufflink_minor = CUFFLINK_MINOR_VERSION;
    int cufflink_subminor = CUFFLINK_SUBMINOR_VERSION;

    std::cout << "Thrust v" << thrust_major << "." << thrust_minor << std::endl;
    std::cout << "Cusp   v" << cusp_major   << "." << cusp_minor   << "."<<cusp_subminor<< std::endl;
    std::cout << "Cufflink   v" << cufflink_major   << "." << cufflink_minor   << "."<<cufflink_subminor<< std::endl;

    return 0;
}
