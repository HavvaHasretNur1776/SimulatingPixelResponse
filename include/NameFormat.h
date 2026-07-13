/*****************************************************************************\
* (c) Copyright 2000-2023 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

#include <sstream>
#include <string>

namespace MPDigitisation {
  /** @brief convenience wrapper around stringstream for formatting strings
   *
   * This class is meant to help in constructing formatted strings for things
   * like histogram names and titles. A second goal is to be lightweight in
   * terms of dependencies and what users need to learn (arguments to operator()
   * are the same that are accepted for std::cout's operator<<).
   *
   * Here's an example:
   * @code
   * std::vector<TH1D> histograms;
   * for (int iLayer = 0; 6 != iLayer; ++i) {
   *     const std::string name = formatter{}("nHitsInLayer", iLayer);
   *     const std::string title = formatter{}("nHits in Layer ", iLayer,
   *         ";nHits;number of events");
   *     histograms.emplace(name.c_str(), title.c_str(), 25, -0.5, 24.5);
   * }
   * @endcode
   *
   * @author Lucas Foreman <lucas.foreman@postgrad.manchester.ac.uk>
   * @author Manuel Schiller <Manuel.Schiller@glasgow.ac.uk>
   * @date 2023-08-2
   */
  class formatter {
  private:
    /// internal stringstream that does all the hard work
    std::stringstream m_str;

  public:
    /// append to the string that's being constructed
    template <class... ARGS>
    formatter& operator()( ARGS&&... args ) {
      const auto nop = []( ... ) {};
      nop( ( m_str << std::forward<ARGS>( args ), 0 )... );
      return *this;
    }

    /// return constructed string
    std::string str() const { return m_str.str(); }
    /// return constructed string
    operator std::string() const { return str(); }
  };

} // namespace MPDigitisation
