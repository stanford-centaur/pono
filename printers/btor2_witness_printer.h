#include <map>
#include <vector>

#include "smt-switch/smt.h"

#include "utils/logger.h"

namespace cosa {

void print_btor_vals_at_time(const smt::TermVec & vec,
                             const smt::UnorderedTermMap & valmap,
                             unsigned int time)
{
  smt::SortKind sk;
  std::string val;
  std::string elem;
  smt::TermVec store_children(3);
  for (size_t i = 0, size = vec.size(); i < size; ++i)
  {
    sk = vec[i]->get_sort()->get_sort_kind();
    if (sk == smt::BV)
    {
      // TODO: this makes assumptions on format of value from boolector
      //       to support other solvers, we need to be more general
      // Remove the #b prefix
      val = valmap.at(vec[i])->to_string();
      val = val.substr(2, val.length() - 2);
      logger.log(0, "{} {} {}@{}", i, val, vec[i], time);
    }
    else if (sk == smt::ARRAY)
    {
      smt::Term tmp = valmap.at(vec[i]);
      while (tmp->get_op() == smt::Store)
      {
        int num = 0;
        for (auto c : tmp)
        {
          store_children[num] = c;
          num++;
        }

        // TODO: this makes assumptions on format of value from boolector
        //       to support other solvers, we need to be more general
        // Remove the #b prefix
        val = store_children[1]->to_string();
        val = val.substr(2, val.length() - 2);
        elem = store_children[2]->to_string();
        elem = elem.substr(2, elem.length() - 2);

        logger.log(0, "{} [{}] {} {}@{}", i, val, elem, vec[i], time);
        tmp = store_children[0];
      }
    }
    else
    {
      throw CosaException("Unhandled sort kind: " + ::smt::to_string(sk));
    }
  }
}

void print_btor_vals_at_time(const std::map<uint64_t, smt::Term> m,
                             const smt::UnorderedTermMap & valmap,
                             unsigned int time)
{
  smt::SortKind sk;
  std::string val;
  std::string elem;
  smt::TermVec store_children(3);
  for (auto entry : m)
  {
    sk = entry.second->get_sort()->get_sort_kind();
    if (sk == smt::BV)
    {
      // TODO: this makes assumptions on format of value from boolector
      //       to support other solvers, we need to be more general
      // Remove the #b prefix
      val = valmap.at(entry.second)->to_string();
      val = val.substr(2, val.length() - 2);
      logger.log(0, "{} {} {}@{}", entry.first, val, entry.second, time);
    }
    else if (sk == smt::ARRAY)
    {
      smt::Term tmp = valmap.at(entry.second);
      while (tmp->get_op() == smt::Store)
      {
        int num = 0;
        for (auto c : tmp)
        {
          store_children[num] = c;
          num++;
        }

        // TODO: this makes assumptions on format of value from boolector
        //       to support other solvers, we need to be more general
        // Remove the #b prefix
        val = store_children[1]->to_string();
        val = val.substr(2, val.length() - 2);
        elem = store_children[2]->to_string();
        elem = elem.substr(2, elem.length() - 2);

        logger.log(
            0, "{} [{}] {} {}@{}", entry.first, val, elem, entry.second, time);
        tmp = store_children[0];
      }
    }
    else
    {
      throw CosaException("Unhandled sort kind: " + ::smt::to_string(sk));
    }
  }
}

void print_witness_btor(const BTOR2Encoder & btor_enc,
                        const std::vector<smt::UnorderedTermMap> & cex)
{
  const smt::TermVec inputs = btor_enc.inputsvec();
  const smt::TermVec states = btor_enc.statesvec();
  const std::map<uint64_t, smt::Term> no_next_states =
      btor_enc.no_next_states();
  bool has_states_without_next = no_next_states.size();

  logger.log(0, "#0");
  print_btor_vals_at_time(states, cex.at(0), 0);

  for (size_t k = 0, cex_size = cex.size(); k < cex_size; ++k)
  {
    // states without next
    if (k && has_states_without_next)
    {
      logger.log(0, "#{}", k);
      print_btor_vals_at_time(no_next_states, cex.at(k), k);
    }

    // inputs
    if (k + 1 < cex_size)
    {
      logger.log(0, "@{}", k);
      print_btor_vals_at_time(inputs, cex.at(k), k);
    }
  }

  logger.log(0, ".");
}

}  // namespace cosa
