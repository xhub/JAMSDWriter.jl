using BinDeps

@BinDeps.setup

libjamsd = library_dependency("libjamsd", aliases=["jamsd"])

# Download binaries from hosted location
bin_prefix = "https://nullptr.fr/lib"

# TODO with latest Julia
if VERSION < v"0.7"
    iswin   = is_windows()
    islinux = is_linux()
    isapple = is_apple()
else
    iswin   = Sys.iswindows()
    islinux = Sys.islinux()
    isapple = Sys.isapple()
end

# TODO remove all is(win|linux|apple) and just have a switch on sys.ARCH?

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

if isapple
    if Sys.ARCH == :x86_64
        provides(Binaries, URI("$bin_prefix-macosx/libjamsd.tar.xz"), libjamsd, os = :Darwin)
    end
end

@BinDeps.install Dict(:libjamsd => :libjamsd)
