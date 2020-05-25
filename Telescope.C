/*
 * Author and Copyright
 * Johannes Gajdosik, 2020
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include "Telescope.H"

#include <map>

typedef std::map<std::string,Telescope::CreationFunction> CreationFunctionMap;

static CreationFunctionMap &GetCreationFunctionMap(void) {
  static CreationFunctionMap creation_function_map;
  return creation_function_map;
}

Telescope::Ptr
Telescope::Create(const std::string &type_and_args,
                  OpenedClosedFunction &&opened_closed,
                  PositionFunction &&announce_position,
                  boost::asio::io_context &io_context) {
  if (type_and_args.empty()) return Ptr();
  std::string::size_type i = type_and_args.find(':');
  CreationFunctionMap &m(GetCreationFunctionMap());
  const auto it(m.find(type_and_args.substr(0,i)));
  if (it == m.end()) return Ptr();
  if (i == std::string::npos) {i = type_and_args.length();}
  else {i++;} // ignore ':'
  Ptr rval(it->second(type_and_args.substr(i),
                      std::move(opened_closed),
                      std::move(announce_position),
                      io_context));
  if (rval->initializationOk()) return rval;
  return Ptr();
}

bool Telescope::RegisterCreationFunction(const std::string &type,
                                         Telescope::CreationFunction &&f) {
  return GetCreationFunctionMap()
           .insert(CreationFunctionMap::value_type(type,std::move(f))).second;
}
