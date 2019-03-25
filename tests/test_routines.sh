
#-----------------------------------------------------------------------------
# test_result()
#-----------------------------------------------------------------------------
# check's the status of the last command and prints results.
# program exits with 'false' status if the last command failed.

test_result()
{
status=$?
if [ "$status" -eq "0" ]; then
    echo "Success"
else
    echo "Failed: $status"
    exit 1
fi
}
