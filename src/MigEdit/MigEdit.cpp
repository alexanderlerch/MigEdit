
// standard headers

// project headers
#include "MigEditConfig.h"

#include "ErrorDef.h"

#include "MigEdit.h"

static const char*  kCMigEditBuildDate             = __DATE__;


CMigEdit::CMigEdit ()
{
    // this never hurts
    this->resetInstance ();
}


CMigEdit::~CMigEdit ()
{
    this->resetInstance ();
}


Error_t CMigEdit::initInstance()
{
    // allocate memory

    // initialize variables and buffers

    return kNoError;
}

Error_t CMigEdit::resetInstance ()
{
    // reset buffers and variables to default values

    return kNoError;
}
