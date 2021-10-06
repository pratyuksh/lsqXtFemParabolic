#include "base_observer.hpp"

#include <fstream>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;


//! Uses config to set class variables.
mymfem::BaseObserver :: BaseObserver (const nlohmann::json& config, int lx)
{
    m_boolVisualize = false;
    if (config.contains("visualization")) {
        m_boolVisualize = config["visualization"];
    }

    m_boolDumpOut = false;
    if (config.contains("dump_output")) {
        m_boolDumpOut = config["dump_output"];
    }

    m_meshNameSuffix = "_lx"+std::to_string(lx);
}

//! Uses GLVis, refer to examples in MFEM documentation.
void mymfem::BaseObserver :: visualize (GridFunction& u) const
{
    socketstream sout;
    if (m_boolVisualize)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        if (!sout)
        {
            std::cout << "Unable to connect to GLVis server at "
                      << vishost << ':' << visport << std::endl;
            m_boolVisualize = false;
            std::cout << "GLVis visualization disabled.\n";
        }
        else
        {
            sout.precision(m_precision);
            sout << "solution\n" << *(u.FESpace()->GetMesh()) << u;
            sout << "pause\n";
            sout << std::flush;
        }
    }
}

//! Writes the mesh to a file.
void mymfem::BaseObserver :: dumpMesh (std::shared_ptr<Mesh>& mesh) const
{
    if (m_boolDumpOut)
    {
        std::string meshName = m_outputDir+"mesh"+m_meshNameSuffix;
        std::cout << meshName << std::endl;

        std::ofstream meshOfs(meshName.c_str());
        meshOfs.precision(m_precision);
        mesh->Print(meshOfs);
        meshOfs.close();
    }
}


// End of file
