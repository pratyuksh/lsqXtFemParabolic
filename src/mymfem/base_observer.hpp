#ifndef BASE_OBSERVER_HPP
#define BASE_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../core/config.hpp"


namespace mymfem {

/**
 * Provides functions to visualize solutions and write output to files.
 */
class BaseObserver
{
public:
    /**
     * @brief Default constructor.
     */
    BaseObserver () {}

    /**
     * @brief Custom constructor.
     * @param config JSON configuration file.
     * @param lx spatial mesh discretisation level.
     */
    BaseObserver (const nlohmann::json&, int);

    /**
     * @brief  Calls the visualize function.
     * @param u scalar grid function, passed by reference.
     */
    void operator () (GridFunction& u) const {
        visualize(u);
    }

    /**
     * @brief Calls the visualize function.
     * @param u scalar grid function, passed as shared pointer.
     */
    void operator() (std::shared_ptr<GridFunction>& u) const {
        visualize(u);
    }

    /**
     * @brief Visualizes the solution.
     * @param u scalar grid function, passed as shared pointer.
     */
    void visualize (std::shared_ptr<GridFunction>& u) const{
        visualize(*u);
    }

    /**
     * @brief Visualizes the solution.
     * @param u scalar grid function, passed by reference.
     */
    void visualize (GridFunction&) const;

    /**
     * @brief Writes mesh to file.
     * @param mesh Mesh passed as a shared pointer.
     */
    void dumpMesh (std::shared_ptr<Mesh>&) const;
    
protected:
    //! Numerical precision of output data.
    int m_precision = 8;

    //! Boolean flag to control visualize.
    mutable bool m_boolVisualize;
    
    //! Boolean flag to control output write.
    bool m_boolDumpOut;
    //! Name of output directory.
    mutable std::string m_outputDir;
    //! Prefix of output mesh file name.
    std::string m_meshNamePrefix;
    //! Suffix of output mesh file name.
    std::string m_meshNameSuffix;
};

}

#endif // BASE_OBSERVER_HPP
