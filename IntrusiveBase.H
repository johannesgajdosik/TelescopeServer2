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

#ifndef INTRUSIVE_BASE_H_
#define INTRUSIVE_BASE_H_

#include <boost/atomic.hpp>

class IntrusiveBase {
protected:
  IntrusiveBase(void) : refcount(0) {}
  virtual ~IntrusiveBase(void) {}
private:
  mutable boost::atomic<int> refcount;
  friend void intrusive_ptr_add_ref(const IntrusiveBase *x) {
    x->refcount.fetch_add(1,boost::memory_order_relaxed);
  }
  friend void intrusive_ptr_release(const IntrusiveBase *x) {
    if (x->refcount.fetch_sub(1,boost::memory_order_release) == 1) {
      boost::atomic_thread_fence(boost::memory_order_acquire);
      delete x;
    }
  }
  IntrusiveBase(const IntrusiveBase&);
  const IntrusiveBase &operator=(const IntrusiveBase&);
};

#endif
