using BinaryProvider

using BinDeps

@BinDeps.setup

libjamsd = library_dependency("libjamsd", aliases=["jamsd"])

# Download binaries from hosted location
bin_prefix = "https://nullptr.fr/lib"

download_info = Dict(
Linux(:x86_64, :glibc)  => ("$bin_prefix-x86_64-linux-gnu/libjamsd.so", "c073dda4e6282aa6ef6cf2f901a86b771fda07412cfc1be735d093446695c47e"),
Windows(:i686)          => ("$bin_prefix-win32/jamsd.dll", "aba26dcb218273545f4aff706a80da90ae74ab32d6aec17ca91efafa5bcf3b29"),
Windows(:x86_64)        => ("$bin_prefix-win64/jamsd.dll", "bee639c681701e2227b939ca47b8eac95b869bd4b0428f83910f7653de6c0e94"),
#    MacOS(:x86_64)          => ("$bin_prefix/libfoo.x86_64-apple-darwin14.tar.gz", "fcc268772d6f21d65b45fcf3854a3142679b78e53c7673dac26c95d6ccc89a24"),
)

# TODO with latest Julia
if VERSION < v"0.7"
	iswin = is_windows()
    islinux = is_linux()
else
	iswin = Sys.iswindows()
    islinux = Sys.is_linux()
end

if iswin
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-win64/libjamsd.tar.xz"), libjamsd, os = :Windows)
    elseif Sys.ARCH == :i686
        provides(Binaries, URI("$bin_prefix-win32/libjamsd.tar.xz"), libjamsd, os = :Windows)
    end
end

if islinux
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-x86_64-linux-gnu/libjamsd.tar.xz"), libjamsd, os = :Linux)
    end
end

@BinDeps.install Dict(:libjamsd => :libjamsd)
