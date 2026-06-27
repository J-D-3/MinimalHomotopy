// Copyright Ingo Proff 2016-2026.
// https://github.com/CrikeeIP/MinimalHomotopy
// Distributed under the MIT Software License (X11 license).

#pragma once

#include <optional>
#include <string>
#include <utility>

namespace mh {

// A minimal "value-or-error" return type — the modern stand-in for the old
// fplus::result. (std::expected would do the same but needs C++23; this keeps
// the project on C++17, which is CGAL's baseline.)
template <class T>
struct Result {
  std::optional<T> value;
  std::string      error;

  bool        ok() const noexcept { return value.has_value(); }
  const T&    operator*() const { return *value; }

  static Result success(T v) { return Result{std::move(v), {}}; }
  static Result failure(std::string e) { return Result{std::nullopt, std::move(e)}; }
};

} // namespace mh
