#pragma once

/** Enum for communicating result of a refinement step.
 *  Used by provers that use abstraction refinement.
 */
enum RefineResult
{
  REFINE_NONE = 0,  // no refinement necessary (e.g. concrete)
  REFINE_SUCCESS,   // refinement successful
  REFINE_FAIL       // failed to refine
};
