/*	MCM file compressor

  Copyright (C) 2015, Google Inc.
  Authors: Mathieu Chartier

  LICENSE

    This file is part of the MCM file compressor.

    MCM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MCM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MCM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "CM.hpp"

namespace cm {

template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::init() {
  const auto start = clock();
  // Simple model.
  {
    size_t idx = 0;
    simple_profile_ = CMProfile();
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder0);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder3);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder6);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder7);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder8);
    if (kInputs > idx++) simple_profile_.EnableModel(kModelOrder9);
    simple_profile_.SetMatchModelOrder(8);
    simple_profile_.SetMinLZPLen(lzp_enabled_ ? 10 : kMaxMatch + 1);
  }
  // Text model.
  const size_t text_mm_order = 7;
  const size_t text_min_lzp_len = lzp_enabled_ ? 12 : kMaxMatch + 1;
  {
    size_t idx = 0;
    text_profile_ = CMProfile();

    if (false) {
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) text_profile_.EnableModel(kModelBracket);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) text_profile_.EnableModel(kModelInterval);
    if (kInputs > idx++) text_profile_.EnableModel(kModelSpecialChar);
    } else {
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) text_profile_.EnableModel(kModelBracket);
    if (kInputs > idx++) text_profile_.EnableModel(kModelInterval);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder3);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder5);
    }
    if (kInputs > idx++) text_profile_.EnableModel(kModelWord1);
    if (kInputs > idx++) text_profile_.EnableModel(kModelInterval2);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder0);
    if (kInputs > idx++) text_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) text_profile_.EnableModel(kModelSparse2);
    if (kInputs > idx++) text_profile_.EnableModel(kModelSparse3);
    if (kInputs > idx++) text_profile_.EnableModel(kModelSparse4);
    if (kInputs > idx++) text_profile_.EnableModel(kModelSparse34);
    // text_profile_ = CMProfile();
    text_profile_.SetMatchModelOrder(text_mm_order);
    text_profile_.SetMinLZPLen(text_min_lzp_len);
  }
  {
    // Text model for match.
    size_t idx = 0;
    text_match_profile_ = CMProfile();

    if (false) {
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelBracket);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelInterval);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSpecialChar);
    } else {
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelBracket);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder7);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelInterval);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSpecialChar);
    }
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelWord1);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder5);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSparse2);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSparse3);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSparse4);
    if (kInputs > idx++) text_match_profile_.EnableModel(kModelSparse34);
    // text_match_profile_ = CMProfile();
    text_match_profile_.SetMatchModelOrder(text_mm_order);
    text_match_profile_.SetMinLZPLen(text_min_lzp_len);

  }
  // Binary model.
  size_t binary_mm_order = 6;
  {
    size_t idx = 0;
    binary_profile_ = CMProfile();
    if (kInputs > idx++) binary_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelOrder3);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelSparse34);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelSparse2);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelInterval2);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelSparse3);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelSparse4);
    if (kInputs > idx++) binary_profile_.EnableModel(kModelOrder0);
    // if (kInputs > idx++) binary_profile_.EnableModel(static_cast<ModelType>(opts_[0]));
    // binary_profile_ = CMProfile();
    binary_profile_.SetMatchModelOrder(binary_mm_order);
    binary_profile_.SetMinLZPLen(lzp_enabled_ ? 0 : kMaxMatch + 1);
    binary_profile_.SetMissFastPath(25000);
  }
  {
    // Binary model for match.
    size_t idx = 0;
    binary_match_profile_ = CMProfile();
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder1);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder2);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelSparse34);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder4);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelInterval2);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelSparse2);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelSparse3);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelSparse4);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder0);
    if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder9);
    // if (kInputs > idx++) binary_match_profile_.EnableModel(kModelOrder3);
    // if (kInputs > idx++) binary_match_profile_.EnableModel(static_cast<ModelType>(opts_[0]));
    // binary_match_profile_ = CMProfile();
    binary_match_profile_.SetMatchModelOrder(binary_mm_order);
  }
  current_interval_map_ = binary_interval_map_;

  table_.build(opts_);

  const size_t extra_mixer_bits = 4 + mem_level_;
  // const size_t extra_mixer_bits = std::min(opt_var_, static_cast<size_t>(4u)) + mem_level_;
  const size_t mixer_bits = kMixerBits;
  const size_t mixer_shift_bits = (kMixerBits - 15);
  mixers_[0].Init(0x100 << extra_mixer_bits, mixer_bits, 25);

  std::cout << std::endl;
  for (auto& m : mixers_) {
    // std::cout << "Mixers " << m.Size() << " RAM=" << m.Size() * sizeof(CMMixer) << " bytes" << std::endl;
  }

  for (auto& c : mixer_text_learn_) c = 9;
  for (auto& c : mixer_binary_learn_) c = 7;
  size_t zero[12] = {};
  if (true) {
    size_t* tl = zero;
    mixer_text_learn_[kModelOrder0] = 24 + tl[11];
    mixer_text_learn_[kModelOrder1] = 24 + tl[9];
    mixer_text_learn_[kModelOrder2] = 24 + tl[1];
    mixer_text_learn_[kModelOrder3] = 11 + tl[6];
    mixer_text_learn_[kModelOrder4] = 12 + tl[0];
    mixer_text_learn_[kModelOrder5] = 9 + tl[7];
    mixer_text_learn_[kModelOrder7] = 9 + tl[3];
    mixer_text_learn_[kModelBracket] = 15 + tl[2];
    mixer_text_learn_[kModelInterval] = 9 + tl[4];
    mixer_text_learn_[kModelInterval2] = 24 + tl[10];
    mixer_text_learn_[kModelWord1] = 20 + tl[5];
    mixer_text_learn_[kModelSpecialChar] = 13 + tl[8];
  }
#if 1
  if (true) {
    size_t tl[] = {4,2,2,1,0,2,2,3,0,1,0,0,};
    // 4,3,4,3,2,3,2,3,1,0,0,0,
    mixer_binary_learn_[kModelOrder0] = 6 + tl[0];
    mixer_binary_learn_[kModelOrder1] = 6 + tl[1];
    mixer_binary_learn_[kModelOrder2] = 6 + tl[2];
    mixer_binary_learn_[kModelOrder3] = 6 + tl[3];
    mixer_binary_learn_[kModelOrder4] = 6 + tl[4];
    mixer_binary_learn_[kModelSparse34] = 6 + tl[5];
    mixer_binary_learn_[kModelSparse2] = 6 + tl[6];
    mixer_binary_learn_[kModelSparse3] = 6 + tl[7];
    mixer_binary_learn_[kModelSparse4] = 6 + tl[8];
    mixer_binary_learn_[kModelInterval2] = 6 + tl[9];
  } else {
    size_t* tl = zero;
    mixer_binary_learn_[kModelOrder0] = 10 + tl[0];
    mixer_binary_learn_[kModelOrder1] = 9 + tl[1];
    mixer_binary_learn_[kModelOrder2] = 10 + tl[2];
    mixer_binary_learn_[kModelOrder3] = 9 + tl[3];
    mixer_binary_learn_[kModelOrder4] = 8 + tl[4];
    mixer_binary_learn_[kModelSparse34] = 9 + tl[5];
    mixer_binary_learn_[kModelSparse2] = 8 + tl[6];
    mixer_binary_learn_[kModelSparse3] = 9 + tl[7];
    mixer_binary_learn_[kModelSparse4] = 7 + tl[8];
    mixer_binary_learn_[kModelInterval2] = 6 + tl[9];
  }
#endif
  for (auto& s : mixer_skip_) s = 0;

  NSStateMap<kShift> sm;
  sm.build(nullptr);

  sse_.init(257 * 256, &table_);
  sse2_.init(257 * 256, &table_);
  sse3_.init(257 * 256, &table_);
  mixer_sse_ctx_ = 0;
  sse_ctx_ = 0;

  hash_mask_ = ((2 * MB) << mem_level_) / sizeof(hash_table_[0]) - 1;
  hash_alloc_size_ = hash_mask_ + kHashStart + (1 << huffman_len_limit);
  hash_storage_.resize(hash_alloc_size_); // Add extra space for ctx.
  hash_table_ = reinterpret_cast<uint8_t*>(hash_storage_.getData()); // Here is where the real hash table starts

  buffer_.Resize((MB / 4) << mem_level_, sizeof(uint32_t));

  // Match model.
  match_model_.resize(buffer_.Size() >> 1);
  match_model_.init(MatchModelType::kMinMatch, 80U);
  fixed_match_probs_.resize(81U * 2);
  int magic_array[100];
  for (size_t i = 1; i < 100; ++i) magic_array[i] = (kMaxValue / 2) / i;
  for (size_t i = 2; i < fixed_match_probs_.size(); ++i) {
    const size_t len = 4 + i / 2;
    auto delta = magic_array[len];
    if ((i & 1) != 0) {
      fixed_match_probs_[i] = table_.st(kMaxValue - 1 - delta);
    } else {
      fixed_match_probs_[i] = table_.st(delta);
    }
  }

  const bool kUseReorder = true;
  uint8_t binary_reorder[] = { 38,2,3,4,5,15,6,23,7,8,9,10,12,13,14,11,17,18,19,16,20,21,24,22,1,25,26,27,28,29,30,31,33,34,35,36,40,37,39,32,42,41,43,44,45,46,47,64,55,49,54,50,48,51,52,53,56,57,58,59,60,61,62,63,84,65,67,66,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,85,86,87,88,89,90,93,94,95,91,92,96,97,98,99,100,101,116,102,103,104,105,106,107,108,109,110,111,112,113,114,115,117,118,119,121,120,122,123,124,125,126,127,128,129,130,131,143,132,133,134,135,136,137,138,139,140,141,142,144,152,145,146,147,148,149,150,151,153,155,154,156,0,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,175,173,174,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,195,192,193,194,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255, };
  uint8_t text_reorder[] = { 7,14,1,12,3,4,11,15,9,16,5,6,18,13,19,30,45,20,21,22,23,17,8,2,26,10,32,43,36,35,42,29,34,24,25,37,31,33,39,38,0,41,28,40,44,58,46,59,92,27,60,61,91,63,95,47,64,124,94,62,93,96,123,125,72,69,65,67,83,68,66,73,82,70,80,76,71,81,77,87,78,74,79,84,75,48,49,50,51,52,53,54,55,56,57,86,88,97,98,99,100,85,101,90,103,104,89,105,107,102,108,109,110,111,106,113,112,114,115,116,119,118,120,121,117,122,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,151,144,145,146,147,148,149,150,152,153,155,156,157,154,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,239,227,228,229,230,231,232,233,234,235,236,237,238,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255, };

  for (int i = 0; i < 256; ++i) {
    // if (opts_) text_reorder[i] = opts_[i];
    text_reorder_[i] = kUseReorder ? text_reorder[i] : i;
    binary_reorder_[i] = kUseReorder ? binary_reorder[i] : i;
  }
  uint8_t binary_mask_map[] = { 15,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,6,12,12,12,12,12,13,12,12,12,12,12,12,12,12,12,12,12,12,12,12,10,12,12,12,12,12,12,12,12,12,12,12,12,12,9,9,9,12,9,9,9,9,9,9,9,12,9,9,9,9,9,9,9,12,9,9,9,12,9,9,9,12,9,9,9,12,7,7,8,12,7,11,7,7,7,14,7,12,7,7,7,12,7,7,7,7,7,7,7,12,7,7,7,12,7,7,7,1,5,5,14,5,5,5,5,5,4,5,3,1,2,5,5,1,5,1,1,1,5,5,5,1,1,1,1,1,1,1,1,1,7,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,1,1,10,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,7,1,1,10,10,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0, };
  uint8_t small_text_mask[] = { 7,7,7,1,4,7,3,7,7,6,7,6,6,6,6,3,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,4,6,5,5,0,2,5,7,2,5,5,7,5,4,3,3,3,3,3,3,3,3,3,3,3,3,5,7,4,1,4,0,2,2,2,2,2,2,2,2,2,2,0,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,5,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,7,7,7,5,7,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, };
  uint8_t text_mask[] = { 15,0,0,2,15,0,8,3,2,12,13,1,3,0,7,9,12,0,0,0,0,0,0,2,0,6,0,0,9,0,0,0,12,7,14,9,7,11,4,11,10,4,9,14,9,8,7,6,5,5,5,5,5,5,5,5,5,5,14,9,2,15,13,4,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,2,4,4,4,4,4,4,5,4,4,3,3,10,1,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, };
  uint8_t text_mask2[] = { 4,2,0,7,2,0,13,0,0,5,4,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,11,2,10,8,5,6,3,9,14,7,7,3,1,5,15,10,0,0,0,0,0,0,0,0,0,0,1,13,13,8,7,7,14,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,15,6,12,14,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,12,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, };
  reorder_.Copy(text_reorder_);
  if (false) {
    std::ofstream of("of.txt");
    for (size_t i = 0; i < 256; ++i) of << (size_t)reorder_.Backward(i) << ",";
    of << std::endl;
    for (size_t i = 0; i < 256; ++i) {
      int c = reorder_.Backward(i);
      of << i << " " << static_cast<int>(small_text_mask[c]) << "," << static_cast<int>(text_mask[c]) << " '";
      of << (isspace(c) ? ' ' : static_cast<char>(c)) << "' ascii=" << (size_t)c << std::endl;
    }
  }
  for (size_t i = 0; i < 256; ++i) {
    // if (opts_) small_text_mask[i] = opts_[i];
    int ri = reorder_[i];
    text_small_interval_map_[ri] = small_text_mask[i];
    text_interval_map_[ri] = text_mask[i];
    text_interval_map2_[ri] = text_mask2[i];
  }
  for (int i = 0; i < 256; ++i) {
    check(text_interval_map_[reorder_[i]] == text_mask[i]);
    check(text_small_interval_map_[reorder_[i]] == small_text_mask[i]);
  }
  word_model_.Init(reorder_);
  bracket_.Init(reorder_);
  special_char_model_.Init(reorder_);

  reorder_.Copy(binary_reorder_);
  for (size_t i = 0; i < 256; ++i) {
    int ri = reorder_[i];
    binary_small_interval_map_[ri] =
      (i < 1) + (i < 32) + (i < 64) + (i < 128) +
      (i < 255) + (i < 142) + (i < 138);
    // binary_interval_map_[ri] = binary_mask_map[i];
    binary_interval_map_[ri] = binary_small_interval_map_[ri] +
      (i < 140) + (i < 137) + (i < 97) + (i < 100) +
      (i < 128) + (i < 48) + (i < 252) + (i < 99);
  }
  // SetActiveReorder(text_reorder);
  current_interval_map_ = binary_interval_map_;
  current_small_interval_map_ = binary_small_interval_map_;
  for (size_t i = 0; i <= WordModel::kMaxLen; ++i) {
    word_model_ctx_map_[i] = (i >= 1) + (i >= 2) + (i >= 3) + (i >= 4) + (i >= 5) + (i >= 6) + (i >= 8);
  }

  // Optimization
  for (uint32_t i = 0; i < kNumStates; ++i) {
    for (uint32_t j = 0; j < 2; ++j) {
      state_trans_[i][j] = sm.getTransition(i, j);
    }
  }

  unsigned short initial_probs[][256] = {
    {1895,1286,725,499,357,303,156,155,154,117,107,117,98,66,125,64,51,107,78,74,66,68,47,61,56,61,77,46,43,59,40,41,28,22,37,42,37,33,25,29,40,42,26,47,64,31,39,0,0,1,19,6,20,1058,391,195,265,194,240,132,107,125,151,113,110,91,90,95,56,105,300,22,831,997,1248,719,1194,159,156,1381,689,581,476,400,403,388,372,360,377,1802,626,740,664,1708,1141,1012,973,780,883,713,1816,1381,1621,1528,1865,2123,2456,2201,2565,2822,3017,2301,1766,1681,1472,1082,983,2585,1504,1909,2058,2844,1611,1349,2973,3084,2293,3283,2350,1689,3093,2502,1759,3351,2638,3395,3450,3430,3552,3374,3536,3560,2203,1412,3112,3591,3673,3588,1939,1529,2819,3655,3643,3731,3764,2350,3943,2640,3962,2619,3166,2244,1949,2579,2873,1683,2512,1876,3197,3712,1678,3099,3020,3308,1671,2608,1843,3487,3465,2304,3384,3577,3689,3671,3691,1861,3809,2346,1243,3790,3868,2764,2330,3795,3850,3864,3903,3933,3963,3818,3720,3908,3899,1950,3964,3924,3954,3960,4091,2509,4089,2512,4087,2783,2073,4084,2656,2455,3104,2222,3683,2815,3304,2268,1759,2878,3295,3253,2094,2254,2267,2303,3201,3013,1860,2471,2396,2311,3345,3731,3705,3709,2179,3580,3350,2332,4009,3996,3989,4032,4007,4023,2937,4008,4095,2048,},
    {2065,1488,826,573,462,381,254,263,197,158,175,57,107,95,95,104,89,69,76,86,83,61,44,64,49,53,63,46,80,29,57,28,55,35,41,33,43,42,37,57,20,35,53,25,11,10,29,16,16,9,27,15,17,1459,370,266,306,333,253,202,152,115,151,212,135,142,148,128,93,102,810,80,1314,2025,2116,846,2617,189,195,1539,775,651,586,526,456,419,400,335,407,2075,710,678,810,1889,1219,1059,891,785,933,859,2125,1325,1680,1445,1761,2054,2635,2366,2499,2835,2996,2167,1536,1676,1342,1198,874,2695,1548,2002,2400,2904,1517,1281,2981,3177,2402,3366,2235,1535,3171,2282,1681,3201,2525,3405,3438,3542,3535,3510,3501,3514,2019,1518,3151,3598,3618,3597,1904,1542,2903,3630,3655,3671,3761,2054,3895,2512,3935,2451,3159,2323,2223,2722,3020,2033,2557,2441,3333,3707,1993,3154,3352,3576,2153,2849,1992,3625,3629,2459,3643,3703,3703,3769,3753,2274,3860,2421,1565,3859,3877,2580,2061,3781,3807,3864,3931,3907,3924,3807,3835,3852,3910,2197,3903,3946,3962,3975,4068,2662,4062,2662,4052,2696,2080,4067,2645,2424,2010,2325,3186,1931,2033,2514,831,2116,2060,2148,1988,1528,1034,938,2016,1837,1916,1512,1536,1553,2036,2841,2827,3000,2444,2571,2151,2078,4067,4067,4063,4079,4077,4075,3493,4081,4095,2048,},
    {1910,1427,670,442,319,253,222,167,183,142,117,119,118,95,82,50,88,92,71,57,53,56,58,52,58,57,32,47,71,37,37,44,42,43,30,25,22,44,16,21,28,64,15,53,27,24,24,12,7,41,28,8,11,1377,414,343,397,329,276,233,200,190,194,230,178,161,157,133,122,110,1006,139,1270,1940,1896,871,2411,215,255,1637,860,576,586,531,573,407,465,353,320,2027,693,759,830,1964,1163,1078,919,923,944,703,2011,1305,1743,1554,1819,2005,2562,2213,2577,2828,2864,2184,1509,1725,1389,1359,1029,2409,1423,2011,2221,2769,1406,1234,2842,3177,2267,3392,2201,1607,3069,2339,1684,3275,2443,3346,3431,3444,3558,3382,3482,3425,1811,1558,3048,3603,3603,3486,1724,1504,2796,3632,3716,3647,3709,2010,3928,2231,3865,2188,3083,2329,2202,2520,2953,2157,2497,2367,3480,3727,1990,3121,3313,3536,2251,2838,2068,3694,3517,2316,3656,3637,3679,3800,3674,2215,3807,2371,1565,3879,3785,2440,2056,3853,3849,3850,3931,3946,3955,3807,3819,3902,3926,2196,3906,3978,3947,3964,4058,2636,4050,2637,4071,2692,2176,4063,2627,2233,1749,2178,2683,1561,1526,2220,947,1973,1801,1902,1652,1434,843,675,1630,1784,1890,1413,1368,1618,1703,2574,2651,2421,2088,2120,1785,2026,4055,4057,4069,4063,4082,4070,3234,4062,4094,2048,},
    {2059,1591,559,339,374,264,145,105,137,137,113,124,89,76,105,69,59,61,41,74,54,46,39,61,36,22,24,57,49,52,41,42,19,37,36,16,21,63,16,50,14,25,33,15,25,23,56,14,41,25,19,20,13,1378,495,348,498,406,255,305,240,237,202,221,187,177,208,164,115,143,1141,178,1250,1984,1770,933,2323,220,259,1482,827,568,666,542,455,466,463,375,399,1998,709,753,811,1770,1225,1182,934,845,873,795,1878,1411,1724,1540,1831,2137,2448,2288,2549,2774,2899,2147,1595,1792,1337,1326,1090,2446,1533,1917,2166,2716,1451,1257,2764,3052,2192,3290,2226,1538,3140,2297,1851,3263,2449,3282,3337,3473,3549,3339,3556,3311,1969,1639,2901,3429,3552,3539,1659,1449,2781,3589,3654,3651,3729,2009,3843,2082,3844,2049,3033,2384,2209,2469,2874,2024,2558,2329,3586,3690,1944,3195,3289,3547,2334,2926,2180,3595,3566,2560,3595,3614,3682,3774,3709,2318,3769,2519,1656,3865,3858,2350,2069,3780,3802,3893,3846,3982,3933,3806,3853,3913,3949,2322,3971,3971,3992,3977,4063,2625,4055,2608,4057,2665,2001,4064,2640,2131,1699,2123,2618,1591,1392,2225,1075,1879,1758,1834,1658,1410,804,688,1600,1852,2015,1685,1551,1629,1527,2435,2387,2210,1995,2023,1816,1711,4062,4048,4050,4084,4063,4078,3311,4087,4095,2048,},
    {2053,1468,621,350,253,230,171,142,158,137,103,145,133,95,88,71,43,80,98,55,56,61,71,32,47,42,61,44,64,58,42,71,53,40,37,50,33,31,40,50,33,41,36,52,19,49,41,22,33,17,14,24,31,1515,550,403,470,433,374,362,248,286,233,323,251,198,194,188,148,206,953,122,1485,2197,1844,900,2585,277,215,1556,779,615,621,474,531,468,423,362,408,1887,732,859,945,1751,1141,1067,944,898,874,827,1897,1434,1687,1547,1970,2091,2394,2273,2407,2650,2800,2175,1743,1792,1483,1415,1007,2307,1796,1966,2150,2664,1658,1451,2632,3041,2183,3301,2170,1498,2999,2249,1843,3304,2498,3368,3352,3476,3499,3377,3559,3182,1964,1888,2852,3462,3489,3505,1819,1221,2735,3547,3643,3659,3692,2179,3828,1987,3847,1935,3153,2245,2185,2242,2730,2101,2765,2617,3576,3760,2085,3019,3170,3437,2529,2931,2149,3559,3553,2645,3608,3620,3684,3778,3673,2334,3798,2607,1520,3834,3817,2530,2147,3749,3782,3862,3801,3905,3928,3835,3936,3912,3981,2355,4004,3956,3958,4019,4046,2497,4009,2571,4038,2574,1910,4039,2602,2191,1783,2207,2882,1689,1505,2423,1110,2271,1854,1908,1848,1341,949,796,1767,1749,2199,1705,1588,1824,1842,2575,2532,2408,2125,2224,1904,1720,4087,4073,4074,4062,4078,4070,3348,4072,4095,2048,},
    {2084,1461,630,443,264,243,221,190,164,167,122,120,134,75,126,88,97,102,145,62,60,80,109,47,97,75,66,89,73,59,76,59,53,68,67,31,53,56,38,11,31,60,54,39,35,51,16,26,33,13,12,21,50,1439,595,447,492,450,354,316,286,253,208,343,276,226,202,141,150,170,756,104,1178,2141,2050,923,2665,287,202,1634,839,631,691,486,518,480,468,414,370,2020,850,832,884,1787,1192,1161,1011,885,902,758,1817,1455,1634,1617,1886,2059,2354,2340,2403,2665,2782,2163,1773,1870,1545,1449,1108,2245,1765,2002,2106,2623,1679,1474,2604,3028,2176,3308,2215,1626,3066,2318,1870,3252,2463,3323,3346,3424,3497,3390,3536,3216,1961,1842,2877,3462,3557,3501,1761,1173,2731,3527,3645,3661,3748,2165,3792,1949,3895,1884,3080,2346,2157,2194,2664,1977,2828,2613,3562,3779,2221,3022,3146,3408,2585,2942,2311,3611,3532,2736,3578,3545,3732,3695,3679,2376,3755,2676,1454,3825,3860,2804,2097,3763,3737,3830,3832,3881,3917,3832,3879,3887,3950,2306,3963,3905,3960,3975,4027,2522,4027,2515,4034,2587,1894,4040,2519,2212,1862,2288,3082,1842,1741,2598,1069,2414,1966,1972,1960,1383,962,993,1922,1843,2264,1467,1434,1684,2075,2831,2781,2743,2412,2407,2090,1955,4084,4085,4092,4079,4082,4088,3658,4083,4095,2048,},
    {2106,1571,503,376,243,143,106,181,107,103,86,81,66,92,79,47,49,68,89,97,49,38,63,60,35,52,65,27,39,46,44,39,37,30,43,42,41,26,25,39,21,26,10,23,25,36,46,24,17,15,11,10,9,1102,414,254,423,342,302,220,245,165,204,250,174,131,166,100,124,152,410,76,1183,1391,1836,839,1641,162,213,1461,791,656,555,443,493,437,439,407,343,1754,744,673,839,1759,1417,1115,1046,983,974,817,1900,1495,1875,1602,1861,2144,2460,2187,2468,2884,2887,2272,1784,1729,1395,1307,1019,2544,1400,1775,2214,2763,1616,1241,2823,3273,2277,3306,2345,1790,3142,2593,2017,3238,2590,3314,3399,3390,3470,3448,3525,3427,2256,1696,3071,3478,3624,3531,1949,1750,2745,3622,3698,3644,3737,2369,3919,2601,3970,2452,3073,2151,2179,2479,2660,1748,2587,2109,3606,3878,2263,2790,3108,3401,2257,2828,2077,3482,3429,2347,3558,3587,3552,3679,3699,2178,3796,2402,1565,3844,3884,2782,2055,3784,3756,3827,3906,3917,3907,3915,3965,3927,4017,2069,3996,4005,4003,3969,4079,2655,4073,2573,4055,2712,2083,4074,2587,2371,2294,2267,3486,2107,2788,2183,1566,2605,2683,2611,1919,1839,1724,1811,2726,2248,1971,2383,2230,2225,2706,3489,3383,3408,2149,3095,2434,2097,4009,3996,3989,4032,4007,4023,2937,4008,4095,2048,},
    {1977,1521,731,469,506,321,288,197,190,139,144,119,169,181,108,135,77,96,63,84,61,57,67,86,33,86,42,70,64,70,74,46,33,31,34,61,29,54,30,26,53,44,28,48,28,58,50,7,26,16,22,28,26,1185,496,372,461,364,305,209,255,243,279,292,202,200,198,154,180,155,551,54,1041,1505,1869,781,1749,289,251,1439,922,725,635,570,604,580,492,414,380,1902,763,803,902,1813,1395,1311,1122,1017,1069,966,1846,1486,1739,1535,1841,2022,2438,2237,2463,2709,2739,2295,1809,1730,1488,1241,1112,2272,1601,1839,2084,2595,1587,1517,2816,3064,2257,3397,2309,1867,3132,2534,2035,3183,2692,3348,3367,3339,3495,3388,3528,3485,2086,1683,3148,3509,3683,3558,2008,1612,2961,3602,3664,3685,3668,2167,3841,2558,3870,2587,3072,2152,2067,2435,2862,2092,2575,1983,3379,3716,2084,3023,3311,3511,2375,3030,2155,3608,3611,2509,3596,3641,3683,3714,3674,2310,3779,2613,1651,3839,3860,2715,2155,3769,3802,3890,3934,3915,3924,3714,3841,3856,3900,2172,3922,3885,3924,3999,4015,2798,4005,2788,4004,2633,2133,4000,2810,2295,1763,2062,3069,1694,2062,2103,1226,2346,2008,2150,1728,1436,1039,1082,1962,1656,1906,2112,2053,1924,2119,3190,3197,3038,1880,2626,2009,1999,4009,3996,3989,4032,4007,4023,2937,4008,4095,2048,},
    {1977,1521,731,469,506,321,288,197,190,139,144,119,169,181,108,135,77,96,63,84,61,57,67,86,33,86,42,70,64,70,74,46,33,31,34,61,29,54,30,26,53,44,28,48,28,58,50,7,26,16,22,28,26,1185,496,372,461,364,305,209,255,243,279,292,202,200,198,154,180,155,551,54,1041,1505,1869,781,1749,289,251,1439,922,725,635,570,604,580,492,414,380,1902,763,803,902,1813,1395,1311,1122,1017,1069,966,1846,1486,1739,1535,1841,2022,2438,2237,2463,2709,2739,2295,1809,1730,1488,1241,1112,2272,1601,1839,2084,2595,1587,1517,2816,3064,2257,3397,2309,1867,3132,2534,2035,3183,2692,3348,3367,3339,3495,3388,3528,3485,2086,1683,3148,3509,3683,3558,2008,1612,2961,3602,3664,3685,3668,2167,3841,2558,3870,2587,3072,2152,2067,2435,2862,2092,2575,1983,3379,3716,2084,3023,3311,3511,2375,3030,2155,3608,3611,2509,3596,3641,3683,3714,3674,2310,3779,2613,1651,3839,3860,2715,2155,3769,3802,3890,3934,3915,3924,3714,3841,3856,3900,2172,3922,3885,3924,3999,4015,2798,4005,2788,4004,2633,2133,4000,2810,2295,1763,2062,3069,1694,2062,2103,1226,2346,2008,2150,1728,1436,1039,1082,1962,1656,1906,2112,2053,1924,2119,3190,3197,3038,1880,2626,2009,1999,4009,3996,3989,4032,4007,4023,2937,4008,4095,2048,},
    {1977,1521,731,469,506,321,288,197,190,139,144,119,169,181,108,135,77,96,63,84,61,57,67,86,33,86,42,70,64,70,74,46,33,31,34,61,29,54,30,26,53,44,28,48,28,58,50,7,26,16,22,28,26,1185,496,372,461,364,305,209,255,243,279,292,202,200,198,154,180,155,551,54,1041,1505,1869,781,1749,289,251,1439,922,725,635,570,604,580,492,414,380,1902,763,803,902,1813,1395,1311,1122,1017,1069,966,1846,1486,1739,1535,1841,2022,2438,2237,2463,2709,2739,2295,1809,1730,1488,1241,1112,2272,1601,1839,2084,2595,1587,1517,2816,3064,2257,3397,2309,1867,3132,2534,2035,3183,2692,3348,3367,3339,3495,3388,3528,3485,2086,1683,3148,3509,3683,3558,2008,1612,2961,3602,3664,3685,3668,2167,3841,2558,3870,2587,3072,2152,2067,2435,2862,2092,2575,1983,3379,3716,2084,3023,3311,3511,2375,3030,2155,3608,3611,2509,3596,3641,3683,3714,3674,2310,3779,2613,1651,3839,3860,2715,2155,3769,3802,3890,3934,3915,3924,3714,3841,3856,3900,2172,3922,3885,3924,3999,4015,2798,4005,2788,4004,2633,2133,4000,2810,2295,1763,2062,3069,1694,2062,2103,1226,2346,2008,2150,1728,1436,1039,1082,1962,1656,1906,2112,2053,1924,2119,3190,3197,3038,1880,2626,2009,1999,4009,3996,3989,4032,4007,4023,2937,4008,4095,2048,},
  };

  for (uint32_t j = 0; j < kProbCtx;++j) {
    for (uint32_t k = 0; k < kNumStates; ++k) {
      int p = initial_probs[std::min(j, 9U)][k];
      probs_[j].SetP(k, p, table_);
      if (j == 0) fast_probs_[k] = table_.st(p);
    }
  }

  /*
  for (size_t i = 0; i < 256; ++i) {
    for (size_t j = 0; j < 256; ++j) {
      fast_mix_[i][j].setP(table_.sq((table_.st(initial_probs[0][i]) + table_.st(initial_probs[1][i])) / 2));
    }
  }*/
  for (size_t i = 0; i < 256 * 256; ++i) {
    fast_mix_[i].setP(2048);
  }
  SetDataProfile(data_profile_);
  last_bytes_ = 0;
  SetUpCtxState();
  // Statistics
  if (kStatistics) {
    for (auto& c : mixer_skip_) c = 0;
    other_count_ = match_count_ = non_match_count_ = 0;
    std::cout << "Setup took: " << clock() - start << std::endl;
    lzp_bit_match_bytes_ = lzp_bit_miss_bytes_ = lzp_miss_bytes_ = normal_bytes_ = 0;
    for (auto& len : match_hits_) len = 0;
    for (auto& len : match_miss_) len = 0;
    miss_len_ = 0;
    for (auto& c : miss_count_) c = 0;
    fast_bytes_ = 0;
  }
}

template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::compress(Stream* in_stream, Stream* out_stream, uint64_t max_count) {
  BufferedStreamWriter<4 * KB> sout(out_stream);
  BufferedStreamReader<4 * KB> sin(in_stream);
  assert(in_stream != nullptr);
  assert(out_stream != nullptr);
  Detector detector(in_stream);
  if (!force_profile_) {
    detector.setOptVar(opt_var_);
    detector.init();
  }
  init();
  ent = Range7();
  if (use_huffman) {
    const clock_t start = clock();
    size_t freqs[256] = { 1 };
    std::cout << "Building huffman tree" << std::endl;
    Huffman::HuffTree* tree = Huffman::Tree<uint32_t>::BuildPackageMerge(freqs, 256, huffman_len_limit);
    tree->PrintRatio(std::cout, "LL");
    Huffman::writeTree(ent, sout, tree, 256, huffman_len_limit);
    huff.build(tree);
    std::cout << "Building huffman tree took: " << clock() - start << " MS" << std::endl;
  }
  for (;max_count > 0; --max_count) {
    uint32_t c;
    if (!force_profile_) {
      Detector::Profile new_profile;
      c = detector.get(new_profile);
      if (new_profile == Detector::kProfileEOF) break;
      auto data_profile = profileForDetectorProfile(new_profile);
      if (data_profile != data_profile_) {
        SetDataProfile(data_profile);
      }
    } else {
      c = sin.get();
      if (c == EOF) break;
    }
    c = reorder_[c];
    dcheck(c != EOF);
    processByte<false>(sout, c);
    update(c);
  }
  ent.flush(sout);

  {
    uint64_t total = 0, less64 = 0;
    for (size_t i = 0; i < 64; ++i) less64 += ctx_count_[i];
    for (size_t i = 0; i < 256; ++i) total += ctx_count_[i];
    if (total > 0) {
      std::cout << std::endl << less64 << "/" << total << " = " << double(less64) / double(total) << std::endl;
    }
  }

  if (kStatistics) {
    if (!kFastStats) {
      std::ofstream fout("probs.txt");
      for (size_t i = 0; i < kInputs; ++i) {
        fout << "{";
        for (uint32_t j = 0; j < 256; ++j) fout << probs_[i].GetP(j) << ",";
        fout << "}," << std::endl;
      }
      // Print average weights so that we find out which contexts are good and which are not.
      for (size_t cur_p = 0; cur_p < static_cast<uint32_t>(kProfileCount); ++cur_p) {
        auto cur_profile = static_cast<DataProfile>(cur_p);
        CMMixer* mixers = mixers_[0].GetMixer();
        std::cout << "Mixer weights for profile " << cur_profile << std::endl;
        for (size_t i = 0; i < 256; ++i) {
          double weights[kInputs + 2] = { 0 };
          double lzp_weights[kInputs + 2] = { 0 };
          size_t count = 0, lzp_count = 0;
          for (size_t j = 0; j < 256; ++j) {
            // Only mixers which have been used at least a few times.
            auto& m = mixers[i * 256 + j];
            if (m.GetLearn() < 30) {
              if (j != 0) {
                for (size_t k = 0; k < m.NumWeights(); ++k) {
                  // weights[k] += double(m.GetWeight(k)) / double(1 << m.shift());
                }
                ++count;
              } else {
                for (size_t k = 0; k < m.NumWeights(); ++k) {
                  // lzp_weights[k] += double(m.GetWeight(k)) / double(1 << m.shift());
                }
                ++lzp_count;
              }
            }
          }
          if (count != 0) {
            std::cout << "Weights " << i << ":";
            for (auto& w : weights) std::cout << w / double(count) << " ";
            std::cout << std::endl;
          }
          if (lzp_count) {
            std::cout << "LZP " << i << ":";
            for (auto& w : lzp_weights) std::cout << w << " ";
            std::cout << std::endl;
          }
        }
      }
      // State count.
      size_t z = 0, nz = 0;
      for (size_t i = 0;i <= hash_mask_;++i) {
        ++(hash_table_[i] != 0 ? nz : z);
      }
      std::cout << "zero=" << z << " nonzero=" << nz << std::endl;
    }
    if (!force_profile_) {
      detector.dumpInfo();
    }
    std::cout << "CMed bytes=" << formatNumber((mixer_skip_[0] + mixer_skip_[1]) / 8)
      << " mix skip=" << formatNumber(mixer_skip_[0])
      << " mix nonskip=" << formatNumber(mixer_skip_[1]) << std::endl;
    std::cout << "match=" << formatNumber(match_count_)
      << " matchfail=" << formatNumber(non_match_count_)
      << " nonmatch=" << formatNumber(other_count_) << std::endl;
    if (!kFastStats) {
      if (false)
        for (size_t i = 0; i < kMaxMatch; ++i) {
          const size_t t = match_hits_[i] + match_miss_[i];
          if (t != 0) {
            std::cout << i << ":" << match_hits_[i] << "/" << match_miss_[i] << " = " << static_cast<double>(match_hits_[i]) / t << std::endl;
          }
        }
      uint64_t miss_tot;
      for (size_t i = 0; i < kMaxMiss;++i) {
        if (miss_count_[i] != 0) std::cout << "Misses " << i << " " << miss_count_[i] + miss_tot << std::endl;
        miss_tot += miss_count_[i];
      }
    }
    std::cout << "Fast bytes " << fast_bytes_ << std::endl;
    if (cur_profile_.MinLZPLen() > 0) {
      std::cout << "lzp_bit_size=" << formatNumber(lzp_bit_match_bytes_)
        << " lzp_bit_miss_bytes=" << formatNumber(lzp_bit_miss_bytes_)
        << " lzp_miss_bytes=" << formatNumber(lzp_miss_bytes_)
        << " lzp_normal_bytes=" << formatNumber(normal_bytes_) << std::endl;
    }
  }
}

template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::decompress(Stream* in_stream, Stream* out_stream, uint64_t max_count) {
  BufferedStreamReader<4 * KB> sin(in_stream);
  BufferedStreamWriter<4 * KB> sout(out_stream);
  Detector detector(out_stream);
  if (!force_profile_) {
    detector.setOptVar(opt_var_);
    detector.init();
  }
  init();
  ent.initDecoder(sin);
  if (use_huffman) {
    // auto* tree = Huffman::readTree(ent, sin, 256, huffman_len_limit);
    // huff.build(tree);
    // delete tree;
  }
  for (; max_count > 0; --max_count) {
    if (!force_profile_) {
      auto new_profile = detector.detect();
      if (new_profile == Detector::kProfileEOF) {
        break;
      }
      auto cm_profile = profileForDetectorProfile(new_profile);
      if (cm_profile != data_profile_) {
        SetDataProfile(cm_profile);
      }
    }
    size_t c = processByte<true>(sin);
    update(c);
    if (force_profile_) {
      sout.put(reorder_.Backward(c));
    } else {
      detector.put(reorder_.Backward(c));
    }
  }
  if (!force_profile_) {
    detector.flush();
  }
  sout.flush();
  size_t remain = sin.remain();
  if (remain > 0) {
    // Go back all the characters we didn't actually read.
    in_stream->seek(in_stream->tell() - remain);
  }
}

template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline CM<kInputs, kUseSSE, HistoryType>::CM(
  const FrequencyCounter<256>& freq,
  uint32_t mem_level,
  bool lzp_enabled,
  Detector::Profile profile)
  : mem_level_(mem_level)
  , data_profile_(profileForDetectorProfile(profile)) {
  force_profile_ = profile != Detector::kProfileDetect;
  lzp_enabled_ = lzp_enabled;
  opts_ = dummy_opts;
  frequencies_ = freq;
}

// Context map for each context.
template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::SetStates(const uint32_t* remap) {
  bool reached[256] = {};
  size_t count = 0;
  for (size_t bits = 0; bits < 255; ++bits) {
    const auto idx = remap[bits];
    for (size_t bit = 0; bit < 2; ++bit) {
      const size_t next_bits = bits * 2 + bit + 1;
      ctx_state_.SetBits(idx, bits);
      if (next_bits < 256) {
        check(reached[next_bits] == false);
        reached[next_bits] = true;
        ++count;
        ctx_state_.SetNext(idx, bit, remap[next_bits]);
      } else {
        ctx_state_.SetNext(idx, bit, next_bits ^ 0x100);
      }
    }
    size_t top_bit = bits + 1;
    while ((top_bit & (top_bit - 1)) != 0) {
      --top_bit;
    }
    ctx_state_.SetBits(idx, (bits + 1) ^ top_bit);
  }
}
  
template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::SetUpCtxState() {
  if (false) {
    OptimalCtxState();
    return;
  }
  uint32_t bits[256] = {256};
  uint32_t ctx_map[256] = {};
  bits[0] = 0;
  for (size_t i = 0; i < 256; ++i) {
    const uint32_t cur_bits = bits[i];
    check(cur_bits != 256);
    for (size_t bit = 0; bit < 2; ++bit) {
      const size_t next = NextNibbleLeaf(i, bit);
      const auto next_bits = cur_bits * 2 + bit + 1;
      if (next < 256) {
        check(next_bits < 256);
        bits[next] = next_bits;
        ctx_map[next_bits] = next;
      }
    }
  }
  SetStates(ctx_map);
}

template <size_t kInputs, bool kUseSSE, typename HistoryType>
inline void CM<kInputs, kUseSSE, HistoryType>::OptimalCtxState() {
  int64_t cost[256] = {};
  // Fill in corresponding
  // byte = (node * 2 + bit + 2) ^ 256
  // (byte ^ 256) = node * 2 + bit + 2
  // (byte ^ 256) - 2 = node * 2 + bit
  auto* freq = frequencies_.GetFrequencies();
  for (size_t i = 0; i < 256; ++i) {
    cost[(i + 256 - 2) / 2] += freq[i];
  }
  for (size_t i = 0; i < 255; ++i) {
    auto next = i * 2 + 2;
    if (next > 255) {
      check(cost[i] == freq[next ^ 256] + freq[(next + 1) ^ 256]);
    }
  }
  uint32_t ctx_map[256] = {};
  size_t cur_ctx = 0;
  for (size_t i = 0; i < 4; ++i) {
    const size_t kTotalSize = 256 * 256;
    int64_t total[kTotalSize] = {};
    std::fill_n(total, kTotalSize, -1);
    const size_t kRemain = kCacheLineSize;
    const size_t byte_states = OptimalByteStates(cost, total, 0, kRemain);
    check(byte_states == total[256 * 0 + kRemain]);
    std::cerr << "Optimal for cache line " << i << " " << byte_states << "/" << frequencies_.Sum() << std::endl;
    using Pair = std::pair<uint32_t, uint32_t>;
    std::vector<Pair> work;
    work.push_back(Pair(0, kRemain));
    size_t cur_cost = 0;
    while (!work.empty()) {
      auto pair = work.back();
      work.pop_back();
      auto node = pair.first;
      auto remain = pair.second;
      const auto next_a = node * 2 + 1;
      const auto next_b = node * 2 + 2;
      if (cost[node] != -1) {
        cur_cost += cost[node];
        cost[node] = -1;
        --remain;
        ctx_map[node] = cur_ctx++;
      }
      if (next_a < 255) {
        check(next_b < 255);
        size_t best_max = 0, best_index = 256;
        for (size_t j = 0; j <= remain; ++j) {
          auto a = std::max(total[256 * next_a + j], int64_t(0));
          auto b = std::max(total[256 * next_b + (remain - j)], int64_t(0));
          if (a + b >= best_max) {
            best_max = a + b;
            best_index = j;
          }
        }
        if (best_index != 0) {
          work.push_back(Pair(next_a, best_index));
        }
        if (remain - best_index != 0) {
          work.push_back(Pair(next_b, remain - best_index));
        }
      }
    }
    check(byte_states == cur_cost);
  }
  for (size_t i = 1; i < 255; ++i) {
    if (ctx_map[i] == 0) {
      ctx_map[i] = cur_ctx++;
    }
  }
  SetStates(ctx_map);
}

}
