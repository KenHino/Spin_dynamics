module spinchem
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, spinchem!"
  end subroutine say_hello
end module spinchem
