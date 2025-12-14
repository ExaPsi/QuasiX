#!/bin/bash
echo "=== QuasiX Rust Code Quality Check ==="
echo ""

# Check for compilation warnings
echo "1. Checking compilation warnings..."
WARNING_COUNT=$(cargo build --release 2>&1 | grep -c "^warning:")
if [ "$WARNING_COUNT" -eq 0 ]; then
    echo "   ✓ No compilation warnings"
else
    echo "   ✗ Found $WARNING_COUNT compilation warnings"
fi

# Check for clippy warnings
echo ""
echo "2. Checking clippy warnings..."
CLIPPY_COUNT=$(cargo clippy --all-features 2>&1 | grep -c "^warning:")
if [ "$CLIPPY_COUNT" -eq 0 ]; then
    echo "   ✓ No clippy warnings"
else
    echo "   ✗ Found $CLIPPY_COUNT clippy warnings"
fi

# Run tests
echo ""
echo "3. Running tests..."
TEST_OUTPUT=$(cargo test --lib 2>&1)
if echo "$TEST_OUTPUT" | grep -q "test result: ok"; then
    TESTS_PASSED=$(echo "$TEST_OUTPUT" | grep -oE "[0-9]+ passed" | grep -oE "[0-9]+")
    echo "   ✓ All $TESTS_PASSED tests passed"
else
    echo "   ✗ Some tests failed"
fi

# Calculate grade
echo ""
echo "=== FINAL GRADE ==="
TOTAL_ISSUES=$((WARNING_COUNT + CLIPPY_COUNT))
if [ "$TOTAL_ISSUES" -eq 0 ] && echo "$TEST_OUTPUT" | grep -q "test result: ok"; then
    echo "Grade: A (Perfect - Zero warnings, all tests pass)"
elif [ "$TOTAL_ISSUES" -le 5 ]; then
    echo "Grade: B (Good - $TOTAL_ISSUES minor issues)"
elif [ "$TOTAL_ISSUES" -le 20 ]; then
    echo "Grade: C (Acceptable - $TOTAL_ISSUES issues)"
else
    echo "Grade: D (Needs work - $TOTAL_ISSUES issues)"
fi
