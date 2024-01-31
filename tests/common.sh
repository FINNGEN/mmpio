function build_mmpio ()
{
    go build -C ../../src -v -ldflags "-X main.MMPioVersion=$(git describe --tags)" -o mmpio
}
