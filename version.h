#pragma once

#define CUFFLINK_VERSION 100

/*! \def CUFFLINK_MAJOR_VERSION
 *  \brief The preprocessor macro \p CUFFLINK_MAJOR_VERSION encodes the
 *         major version number of the Cufflink library.
 */
#define CUFFLINK_MAJOR_VERSION     (CUFFLINK_VERSION / 100000)

/*! \def CUFFLINK_MINOR_VERSION
 *  \brief The preprocessor macro \p CUFFLINK_MINOR_VERSION encodes the
 *         minor version number of the CUFFLINK library.
 */
#define CUFFLINK_MINOR_VERSION     (CUFFLINK_VERSION / 100 % 1000)

/*! \def CUFFLINK_SUBMINOR_VERSION
 *  \brief The preprocessor macro \p CUFFLINK_SUBMINOR_VERSION encodes the
 *         sub-minor version number of the CUFFLINK library.
 */
#define CUFFLINK_SUBMINOR_VERSION  (CUFFLINK_VERSION % 100)


