#if !defined(__CommandLineOptions_hdr__)
#define __CommandLineOptions_hdr__

#include <string>
#include "ErrorDef.h"

class CCommandLineOptions
{
public:

    class COption
    {
    public:
        enum OptionType_t
        {
            kBool,
            kString,
            kInt,
            kFloat,

            kNumOptionTypes
        };

        COption(int iIdx, std::string sClArg, std::string sDesc, OptionType_t eType) : 
            m_iIdentifier(iIdx),
            m_sSwitch(sClArg),
            m_sDesc(sDesc),
            m_eType(eType){};
        virtual ~COption(){};

        int getIdentifier()
        {
            return m_iIdentifier;
        }
        std::string getClArg()
        {
            return m_sSwitch;
        }
        std::string getDesc()
        {
            return m_sDesc;
        }
        OptionType_t getType()
        {
            return m_eType;
        }
    private:
        int             m_iIdentifier;
        std::string     m_sSwitch;
        std::string     m_sDesc;
        OptionType_t    m_eType;
    };

    /*! creates a new CommandLineOptions instance
    \param CCommandLineOptions * & pCCommandLineOptions: pointer to the new instance
    \return Error_t
    */
    static Error_t createInstance (CCommandLineOptions*& pCCommandLineOptions);
    
    /*! destroys an CommandLineOptions instance
    \param CCommandLineOptions * & pCCommandLineOptions: pointer to the instance to be destroyed
    \return Error_t
    */
    static Error_t destroyInstance (CCommandLineOptions*& pCCommandLineOptions);
    
    /*! initializes an CommandLineOptions instance
    \param argc number of cl arguments
    \param argv arguments
    \return Error_t
    */
    Error_t initInstance (const COption myOptions[], int iNumOfOptions);
    
    /*! resets an CommandLineOptions instance
    \return Error_t
    */
    Error_t resetInstance ();

    /*! parses the command line arguments
    \param int argc
    \param char * argv[]
    \return Error_t
    */
    Error_t process(int argc, char* argv[]);

    /*! returns if the option has been set on the command line
    \param int iIdentifier: integer ID as given in COptions class
    \return bool
    */
    bool isOptionSet(int iIdentifier) const;
    /*! returns if the option has been set on the command line
    \param std::string argv: command line switch to check
    \return bool
    */
    bool isOptionSet( std::string argv ) const;

    /*! returns the option value as string
    \param int iIdentifier: integer ID as given in COptions class
    \param std::string & sResult
    \return Error_t
    */
    Error_t getOption(int iIdentifier, std::string& sResult) const;
    /*! returns the option value as string
    \param std::string argv: command line switch to check
    \param std::string & sResult
    \return Error_t
    */
    Error_t getOption(std::string argv, std::string& sResult) const;

protected:
    CCommandLineOptions ();
    virtual ~CCommandLineOptions () {};
private:
    class CClArg 
    {
    public:
        CClArg():
            m_bIsSet(false)
        {
        };
        virtual ~CClArg(){};

        void setCurrOption(std::string sValue = "")
        {
            m_sOption   = sValue;
            m_bIsSet    = true;
        }
        std::string getArg() const
        {
            if (isSet())
                return m_sOption;
            else
            {
                return "";
            }
        }
        void reset()
        {
            m_bIsSet = false;
        }
        bool isSet() const
        {
            return m_bIsSet;
        }

        bool m_bIsSet;
        std::string m_sOption;
    };

    int getBufferIdx( std::string argv ) const;
    int getBufferIdx( int iIdentifier) const;
    
    bool m_bIsInitialized;

    COption **m_ppCOptions;
    int     m_iNumOfOptions;
    CClArg  **m_ppCClArgs;
};

#endif // #if !defined(__CommandLineOptions_hdr__)



